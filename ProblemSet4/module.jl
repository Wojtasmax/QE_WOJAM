module Consumption_Savings_Model

using QuantEcon, Interpolations, LinearAlgebra, Parameters, Printf, Roots, Optim, Statistics, Printf, Random

export Consumption_Savings_Model

@with_kw struct Consumption_Savings
    #Utiity function
    u = γ == 1 ? (c -> log(c)) : (c -> (c^(1 - γ) - 1) / (1 - γ))
    u_prime = γ == 1 ? (c -> 1/c) : (c -> c^(-γ))
    u_prime_inv = γ == 1 ? (y -> 1/y) : (y -> y^(-1/γ))

    γ = 2.0 
    β = 0.99 
    R = 1.010 

    #AR process for income
    ρ_z = 0.90 #income persistance
    σ_ϵ = 0.2*sqrt(1 - ρ_z^2) 
    N_z = 5  #number of wage levels
    mc_z = rouwenhorst(N_z, ρ_z, σ_ϵ, 0.0)
    λ_z = stationary_distributions(mc_z)[1]
    P_z = mc_z.p
    z_vec = exp.(mc_z.state_values) / sum(exp.(mc_z.state_values) .* λ_z) # normalize mean to 1
    z_min = minimum(mc_z.state_values) #calculating z_min for the borrowing limit

    #asset grid 
    a_min = -0.6 * z_min / (R - 1) 
    a_max = 500.0
    N_a = 100
    θ = 3
    ω = range(0, 1, length = N_a) 
    a_vec = a_min .+ (a_max - a_min) .* ω.^θ 

    L = sum(z_vec .* λ_z) #expected income
end

#The first method is standard value function iteration (VFI) with linear interpolation. 
#VFI with interpolation (direct maximization)   
function T_interp(v, model)
    @unpack N_a, N_z, a_vec, z_vec, P_z, β, u,R , a_min, a_max = model

    v_new = zeros(N_a, N_z)
    σ_new = zeros(N_a, N_z)

    v_interps = [LinearInterpolation(a_vec, v[:, iz], extrapolation_bc = Line()) for ia in 1:N_z]

    for (iz, z) in enumerate(z_vec)
        # Pre-compute expected value function E[v(k', z') | z] for any k'
        function EV(a_next)
            ev = 0.0
            for iz_next in 1:N_z
                ev += P_z[iz, iz_next] * v_interps[iz_next](a_next)
            end
            return ev
        end
        
        for (ia, a) in  enumerate(a_vec)
            cash = R*a + z 
            a_max_feasible = cash - 1e-10
            
            if a_max_feasible < a_min
                σ_new[ia, iz] = a_min



            function objective(a_next)
                c = R*a + z - a_next
                if c <= 0
                    return -Inf
                else
                    return -(u(c) + β * EV(a_next)) #minimization so negative
                end
            end

            if a_high <= a_low 
                σ_new[ia, iz] = a_low
                v_new[ia, iz] = u(wealth - a_low) + β * EV(a_low)
            else
                result = optimize(objective, a_low, a_high, Brent())
                σ_new[ia, iz] = Optim.minimizer(result)
                v_new[ia, iz] = Optim.minimum(result) 
            end
        end
    end

    return v_new, σ_new
end

function vfi_interp(model; maxiter = 1000, tol = 1e-7, v_init = nothing)
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
            



end #ending module

