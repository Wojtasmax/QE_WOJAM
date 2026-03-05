module Engine

using QuantEcon, Statistics, Parameters, Interpolations, Optim, NumericalIntegration, LinearAlgebra
export ProjectParams, InitializeModel, Solve_Model, VFI, GENERATE_K_GRID, ADJ_COST, IRR, OPTIMAL_H, OPERATING_PROFIT, get_transition_matrix_young, stationary_distribution, moments, objective_function_SMM, new_objective_function_SMM

@with_kw struct ProjectParams
    #=====Parametry globalne=====#
    alpha = 0.3
    v = 0.6
    r = 0.04
    delta = 0.08
    w = 1.0
    beta = 1.0 / (1.0 + 0.04)

    #=====Parametry procesu Markowa=====#
    rho::Float64 = 0.9
    sigma_eps::Float64 = 0.12
    N_z::Int = 7  #to jest z czapy - zmienić na 3

    #=====Wyniki procesu Markowa=====#
    z_vec::Vector{Float64} = Float64[]
    P_z::Matrix{Float64} = zeros(1,1)
    lambda_z::Vector{Float64} = Float64[]


    #===Adjustment Cost====#
    #losowe wartosci    
    gamma = 0.05          
    F = 0.01          
    ps = 0.80         

    #======Grid=====#
    # jak bedzie mulic to zmienic 
    #jak beda bledy numeryczne to zmienic k_min
    k_max = 500
    k_min = 1e-4
    N_A = 500   
    omega = range(0, 1, length=N_A)
end


function InitializeModel(; alpha=0.3, v=0.6, r=0.04, delta=0.08, w=1.0, rho=0.9, sigma_eps=0.12, N_z=7, gamma=0.05,
    F=0.01, ps=0.80, k_max=500, k_min=1e-4, N_A=500)
    beta = 1.0 / (1.0 + r)
    omega = range(0, 1, length=N_A)  
    
    z_tilde = exp(-sigma_eps^2 / (2 * (1 - rho^2)))
    mu_logz = log(z_tilde)
    
    mc_z = rouwenhorst(N_z, rho, sigma_eps, mu_logz)
    P_z = mc_z.p
    lambda_z = stationary_distributions(mc_z)[1]
    
    z_raw = exp.(mc_z.state_values)
    z_vec = z_raw ./ sum(z_raw .* lambda_z)
    
    return ProjectParams(
        alpha=alpha, v=v, r=r, delta=delta, w=w, beta=beta, 
        rho=rho, sigma_eps=sigma_eps, N_z=N_z, 
        z_vec=z_vec, P_z=P_z, lambda_z=lambda_z,
        gamma=gamma, F=F, ps=ps, k_min=k_min, k_max=k_max, N_A=N_A, omega=omega
    )
end

function ADJ_COST(model::ProjectParams, i, k, toll_level=0.0005)
    @unpack F, gamma = model
    if  (-toll_level<=i<=toll_level)
        return 0.0
    else
        return ((gamma/2)*(i/k)^2)*k + F*k
    end
end


function IRR(model::ProjectParams, i)
    @unpack ps = model
    if i >= 0
        return 1.0
    else
        return ps
    end
end


function OPTIMAL_H(model::ProjectParams, k, z)
    @unpack w, v, alpha = model
    return (w/(v*z*k^alpha))^(1/(v-1))
end


function OPERATING_PROFIT(model::ProjectParams, k, z)  
    @unpack alpha, v, w = model
    h = OPTIMAL_H(model, k, z)
    return z*k^alpha*h^v - w*h
end

function GENERATE_K_GRID(model::ProjectParams, type=:polynomial, theta=5)
    @unpack omega, k_min, k_max, N_A = model

    if type == :polynomial
        return  k_min .+ (k_max - k_min) .* collect(omega).^theta
    elseif type == :exp
        return exp.(range(log(k_min), log(k_max), length=N_A))
    end
end


function VFI(model::ProjectParams, V_old::Matrix{Float64}, k_grid_type=:polynomial, theta=5.0)
    @unpack alpha, beta, v, r, w, delta, P_z, z_vec, N_A, N_z = model
    K_GRID = GENERATE_K_GRID(model, k_grid_type, theta)

    V_new = zeros(N_A, N_z)
    INVESTMENT_POLICY = zeros(N_A, N_z)
    FUTURE_CAPITAL_POLICY = zeros(N_A, N_z)

    V_expect = V_old * P_z'

    Threads.@threads for z_idx in 1:N_z
        z = z_vec[z_idx]
        for k_idx in 1:N_A
            k = K_GRID[k_idx]
            best_value = -Inf
            OP_PROFIT = OPERATING_PROFIT(model, k, z)

            for knxt_idx in 1:N_A
                k_next = K_GRID[knxt_idx]
                i = k_next - (1-delta)*k
                PROFIT_NOW = OP_PROFIT - ADJ_COST(model, i, k) - IRR(model, i)*i
                value = PROFIT_NOW + beta * V_expect[knxt_idx, z_idx]

                if value > best_value
                    V_new[k_idx, z_idx] = value
                    INVESTMENT_POLICY[k_idx, z_idx] = i
                    FUTURE_CAPITAL_POLICY[k_idx, z_idx] = k_next
                    best_value = value
                end
            end
        end
    end
    return V_new, INVESTMENT_POLICY, FUTURE_CAPITAL_POLICY
end
    
function Solve_Model(model::ProjectParams, V_old, min_error, max_iter)
    error = 1.0
    iter = 0
    total_time = @elapsed begin
        while (error > min_error && iter < max_iter)
            out = VFI(model, V_old, :polynomial, 5.0)
            iter += 1
            V_new = out[1]
            error = maximum(abs.(V_new .- V_old))
            V_old = V_new
            if iter % 50 == 0
                println("Iteracja $iter, błąd = $error")
            end
        end
    end
    println("Całkowity czas obliczeń: ", round(total_time, digits=3), " sekund.")
    println("Średni czas na iterację: ", round(total_time/iter, digits=5), " s.")
    return VFI(model, V_old, :polynomial, 5.0)
end


function get_transition_matrix_young(model::ProjectParams, future_k_policy::Matrix{Float64})
    # parameters
    @unpack N_A, N_z, P_z = model
    
    # k grid generation
    K_GRID = GENERATE_K_GRID(model, :polynomial, 5.0)
    
    Q = zeros(N_A * N_z, N_A * N_z)
    
    for iz in 1:N_z
        for ik in 1:N_A
            k_next = future_k_policy[ik, iz]
            
            if k_next <= K_GRID[1]
                ik_low = 1
                ik_high = 1
                weight_low = 1.0
                weight_high = 0.0
            elseif k_next >= K_GRID[end]
                ik_low = N_A
                ik_high = N_A
                weight_low = 1.0
                weight_high = 0.0
            else
                ik_high = searchsortedfirst(K_GRID, k_next)
                ik_low = ik_high - 1
                
                # linear interpolaion weights
                weight_high = (k_next - K_GRID[ik_low]) / (K_GRID[ik_high] - K_GRID[ik_low])
                weight_low = 1.0 - weight_high
            end
            
            # probability distribution for next z state
            for iz_next in 1:N_z
                row = (iz - 1) * N_A + ik
                col_low = (iz_next - 1) * N_A + ik_low
                col_high = (iz_next - 1) * N_A + ik_high
                
                Q[row, col_low] += P_z[iz, iz_next] * weight_low
                Q[row, col_high] += P_z[iz, iz_next] * weight_high
            end
        end
    end
    
    return Q
end

function stationary_distribution(model::ProjectParams, future_k_policy::Matrix{Float64})
    @unpack N_A, N_z = model
    
    Q = get_transition_matrix_young(model, future_k_policy)
    
    # initial guess for the distribution
    λ_vector = ones(N_A * N_z) / (N_A * N_z)
    
    # stationary distribution calculation using power iteration
    for iter in 1:10000
        λ_new = Q' * λ_vector
        
        if maximum(abs.(λ_new - λ_vector)) < 1e-10
            λ_vector = λ_new
            break
        end
        λ_vector = λ_new
    end
    
    # standardize to sum to 1
    λ_vector = λ_vector / sum(λ_vector)
    
    # reshape to (N_A, N_z) for easier interpretation
    μ = zeros(N_A, N_z)
    for iz in 1:N_z
        μ[:, iz] = λ_vector[(iz-1)*N_A+1:iz*N_A]
    end
    
    return μ
end

function moments(μ, result, K_GRID, verbose=false)
    i_grid = result[2]
    i_k = i_grid ./ K_GRID
    average_investment_rate = sum(i_k .* μ) / sum(μ)
    inaction_rate = sum((abs.(i_k) .<= 0.01) .* μ) / sum(μ)
    fraction_with_negative_investment = sum((i_grid .< 0) .* μ) / sum(μ)
    positive_spike_rate = sum((i_k .> 0.2) .* μ) / sum(μ)
    negative_spike_rate = sum((i_k .< -0.2) .* μ) / sum(μ)
    if verbose
        println("\n=== MOMENTY ROZKŁADU ===")
        println("Average Investment Rate: ", round(average_investment_rate, digits=4))
        println("Inaction Rate (|i/k| ≤ 1%): ", round(inaction_rate, digits=4))
        println("Fraction with Negative Investment: ", round(fraction_with_negative_investment, digits=4))
        println("Positive Spike Rate (i/k > 20%): ", round(positive_spike_rate, digits=4))
        println("Negative Spike Rate (i/k < -20%): ", round(negative_spike_rate, digits=4))
    end
    return (average_investment_rate, inaction_rate, fraction_with_negative_investment, 
            positive_spike_rate, negative_spike_rate)
end



function objective_function_SMM(theta::Vector{Float64}, verbose=false)
    gamma_val, F_val, ps_val = theta[1], theta[2], theta[3]
    
    TARGET_MOMENTS = [0.122, 0.081, 0.104, 0.180, 0.014]

    if gamma_val <= 0 || F_val < 0 || ps_val < 0 || ps_val > 1
        return 1e6 
    end

    if verbose
    println("Testuję parametry: γ = $(round(gamma_val, digits=4)), F = $(round(F_val, digits=4)), ps = $(round(ps_val, digits=4))")
    end

    N_z_opt = 3
    rho_opt = 0.9
    sigma_eps_opt = 0.12

    z_tilde = exp(-sigma_eps_opt^2 / (2 * (1 - rho_opt^2)))
    mu_logz = log(z_tilde)
    mc_z = rouwenhorst(N_z_opt, rho_opt, sigma_eps_opt, mu_logz)
    
    P_z_opt = mc_z.p
    lambda_z_opt = stationary_distributions(mc_z)[1]
    z_raw = exp.(mc_z.state_values)
    z_vec_opt = z_raw ./ sum(z_raw .* lambda_z_opt)

    model_opt = ProjectParams(
        gamma=gamma_val, F=F_val, ps=ps_val, 
        N_A=100, N_z=N_z_opt, 
        P_z=P_z_opt, z_vec=z_vec_opt, lambda_z=lambda_z_opt
    )
    # można zmienić 1e-4 by było szybciej/wolniej
    V_init = zeros(model_opt.N_A, model_opt.N_z)
    result_opt = Solve_Model(model_opt, V_init, 1e-4, 2000) 

    μ_opt = stationary_distribution(model_opt, result_opt[3])
    K_GRID_opt = GENERATE_K_GRID(model_opt, :polynomial, 5.0)
    
    moms_tuple = moments(μ_opt, result_opt, K_GRID_opt)
    model_moms = [moms_tuple...] 
    
    distance = sum((model_moms .- TARGET_MOMENTS).^2)
    println("Obecny błąd: ", round(distance, digits=6), "\n")
    
    return distance
end

end #moduł