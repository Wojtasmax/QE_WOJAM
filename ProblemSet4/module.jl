module Consumption_Savings_Model

using QuantEcon, Interpolations, LinearAlgebra, Parameters, Printf, Roots, Optim, Statistics, Random

export Consumption_Savings, T_interp, vfi_interp, create_initial_guess, get_consumption_matrix, euler_residuals, simulate_path


@with_kw struct Consumption_Savings    #Key Parameters
    γ = 2.0
    β = 0.99
    R = 1.010

    #Utility function
    u = γ == 1 ? (c -> log(c)) : (c -> (c^(1 - γ) - 1) / (1 - γ))
    u_prime = γ == 1 ? (c -> 1 / c) : (c -> c^(-γ))
    u_prime_inv = γ == 1 ? (y -> 1 / y) : (y -> y^(-1 / γ))

    #AR process for income
    ρ_z = 0.90 #income persistance
    σ_ϵ = 0.2 * sqrt(1 - ρ_z^2)
    N_z = 3  #number of income levels
    mc_z = rouwenhorst(N_z, ρ_z, σ_ϵ, 0.0)
    λ_z = stationary_distributions(mc_z)[1]
    P_z = mc_z.p
    z_vec = exp.(mc_z.state_values) / sum(exp.(mc_z.state_values) .* λ_z) # normalize mean to 1
    z_min = minimum(z_vec) #calculating z_min for the borrowing limit

    #asset grid 
    a_min = -0.6 * z_min / (R - 1)
    a_max = 500.0
    N_a = 100
    θ = 3
    ω = range(0, 1, length=N_a)
    a_vec = a_min .+ (a_max - a_min) .* ω .^ θ

    L = sum(z_vec .* λ_z) #expected income
end


function create_initial_guess(model)
    # Create a reasonable initial value function based on steady-state consumption
    @unpack N_a, N_z, a_vec, z_vec, u, β, R = model
    
    v_init = zeros(N_a, N_z)
    
    for (iz, z) in enumerate(z_vec)
        for (ia, a) in enumerate(a_vec)

            income = z + R*a  # income at the current state
            c_heuristic = 0.7 * income  # consume 70%, save 30%
            
            if c_heuristic > 0
                # Approximate value as present value of constant consumption
                v_init[ia, iz] = u(c_heuristic) / (1 - β)
            else
                v_init[ia, iz] = -1e10
            end
        end
    end
    
    return v_init
end

#The first method is standard value function iteration (VFI) with linear interpolation. 
#VFI with interpolation (direct maximization)   
function T_interp(v, model)
    @unpack N_a, N_z, a_vec, z_vec, P_z, β, u, R, a_min, a_max = model

    v_new = zeros(N_a, N_z)
    σ_new = zeros(N_a, N_z)

    v_interps = [LinearInterpolation(a_vec, v[:, iz], extrapolation_bc=Line()) for iz in 1:N_z]

    for (iz, z) in enumerate(z_vec)
        # Pre-compute expected value function E[v(k', z') | z] for any k'
        function EV(a_next)
            ev = 0.0
            for iz_next in 1:N_z
                ev += P_z[iz, iz_next] * v_interps[iz_next](a_next)
            end
            return ev
        end

        for (ia, a) in enumerate(a_vec)

            function objective(a_next)
                c = R * a + z - a_next
                if c <= 0
                    return Inf
                else
                    return -(u(c) + β * EV(a_next)) #minimization so negative
                end
            end

            wealth = R * a + z
            a_max_feasible = wealth - 1e-10

            a_low = a_min
            a_high = min(a_max_feasible, a_max)

            if a_high < a_low
                σ_new[ia, iz] = a_low
                v_new[ia, iz] = -Inf

            else
                result = optimize(objective, a_low, a_high, Brent())
                σ_new[ia, iz] = Optim.minimizer(result)
                v_new[ia, iz] = -Optim.minimum(result)
            end
        end
    end

    return v_new, σ_new
end

function vfi_interp(model; maxiter=1000, tol=1e-7, v_init=nothing)
    @unpack N_a, N_z = model

    # Use provided initial guess or create one
    v = isnothing(v_init) ? zeros(N_a, N_z) : copy(v_init)
    σ = zeros(N_a, N_z)
    err = tol + 1.0
    iter = 1

    while err > tol && iter < maxiter
        v_new, σ = T_interp(v, model)
        err = maximum(abs.(v_new - v) ./ (1.0 .+ abs.(v)))
        v = v_new
        iter += 1
    end

    println("VFI with interpolation converged in $iter iterations, error: $err")
    return v, σ, iter, err
end

function get_consumption_matrix(model, policy_a)
    @unpack N_a, N_z, a_vec, z_vec, R = model
    c_mat = zeros(N_a, N_z)
    
    for iz in 1:N_z
        for ia in 1:N_a
            wealth = R * a_vec[ia] + z_vec[iz]
            c_mat[ia, iz] = wealth - policy_a[ia, iz]
        end
    end
    return c_mat
end


function euler_residuals(model, σ; test_grid=nothing)
    # Compute Euler equation residuals to verify solution accuracy.
    # Euler equation: u'(c) = β*R * E[u'(c')]
    # Residual: relative error in consumption implied by Euler equation
    @unpack a_vec, z_vec, P_z, N_z, β, R, γ, u_prime, u_prime_inv = model 
    
    if isnothing(test_grid)
        test_grid = range(model.a_min + 0.01, model.a_max * 0.5, length=500)
    end
    
    n_test = length(test_grid)
    residuals = zeros(n_test, N_z)
    

    σ_interps = [LinearInterpolation(a_vec, σ[:, iz], extrapolation_bc=Line()) for iz in 1:N_z]
    
    for (iz, z) in enumerate(z_vec)
        policy_interp = σ_interps[iz]
        

        for (ia, a) in enumerate(test_grid)
            a_next = policy_interp(a)
            c = R * a + z - a_next
            
            if c > 1e-10
                expected_mu_next = 0.0
                

                for iz_next in 1:N_z
                    z_next = z_vec[iz_next]
                    

                    a_next_next = σ_interps[iz_next](a_next)
                    c_next = R * a_next + z_next - a_next_next
                    
                    if c_next > 1e-10
                         expected_mu_next += P_z[iz, iz_next] * u_prime(c_next)
                    else 
                         # Penalty for non-positive consumption
                         expected_mu_next += P_z[iz, iz_next] * 1e10 
                    end
                end
                

                euler_rhs = β * R * expected_mu_next
                

                c_implied = u_prime_inv(euler_rhs)
                
                residuals[ia, iz] = (c - c_implied) / c
            else
                residuals[ia, iz] = NaN
            end
        end 
    end 
    
    return test_grid, residuals
end 

function simulate_path(model, σ, T::Int = 10000)
    @unpack a_vec, z_vec, P_z, N_z, R = model 

    a_sim = Vector{Float64}(undef, T)
    z_sim = Vector{Float64}(undef, T)
    iz_sim = Vector{Int}(undef, T)
    c_sim = Vector{Float64}(undef, T)

    iz_sim[1] = N_z ÷ 2 + 1
    z_sim[1] = z_vec[iz_sim[1]]
    a_sim[1] = 0

    for t in 1:T-1
        σ_interp_fn = LinearInterpolation(a_vec, σ[:, iz_sim[t]], extrapolation_bc = Line())
        a_sim[t + 1] = σ_interp_fn(a_sim[t])
        c_sim[t] = R*a_sim[t] + z_sim[t] - a_sim[t + 1]

        iz_sim[t + 1] = findfirst(cumsum(P_z[iz_sim[t], :]) .>= rand())
        z_sim[t + 1] = z_vec[iz_sim[t + 1]]
    end
    c_sim[T] = R * a_sim[T] + z_sim[T] - a_sim[T] #consumption in the last period

    return a_sim[], z_sim[], c_sim, iz_sim #not burning the initial 500 here
end


end