module Engine

using QuantEcon, Statistics, Parameters, Interpolations, Optim
export ProjectParams, Solve_Model, VFI, GENERATE_K_GRID, ADJ_COST, IRR, OPTIMAL_H, OPERATING_PROFIT

@with_kw struct ProjectParams
    #=====Parametry globalne=====#
    alpha = 0.3
    v = 0.6
    r = 0.04
    delta = 0.08
    w = 1.0
    beta = 1.0 / (1.0 + r)
    
    #=====Parametry procesu Markowa=====#
    rho::Float64 = 0.9
    sigma_eps::Float64 = 0.12
    N_z::Int = 7  #to jest z czapy

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


function ProjectParams(model::ProjectParams, alpha=0.3, v=0.6, r=0.04, delta=0.08, w=1.0, rho=0.9, sigma_eps=0.12, N_z=7, gamma=0.05,
    F=0.01, ps=0.80, k_max=500, k_min=1e-4, N_A=500)
    @unpack rho, sigma_eps, N_z = model

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
    if  (-toll_level<=i<=toll_level)
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
        return  k_min .+ (k_max - k_min) .* omega.^theta
    elseif type == :exp
        return exp.(range(log(k_min), log(k_max), length=N_A))
    end
end


function VFI(model::ProjectParams, V_old::Matrix{Float64}, k_grid_type=:polynomial, theta=5.0)
    @unpack alpha, beta, v, r, w, delta, P_z, z_vec = model
    K_GRID = GENERATE_K_GRID(model, k_grid_type, theta)

    V_new = copy(V_old)
    INVESTMENT_POLICY = copy(V_old)
    FUTURE_CAPITAL_POLICY = copy(V_old)
    ADJ_COST_POLICY = copy(V_old)
    IRR_COST_POLICY = copy(V_old)
    OPERATING_PROFIT_POLICY = copy(V_old)
    TOTAL_PROFIT_POLICY = copy(V_old)

    for (z_idx, z) in enumerate(z_vec)
        for (k_idx, k) in enumerate(K_GRID)
            best_value = -Inf
            OP_PROFIT = OPERATING_PROFIT(model, k, z)

            for (knxt_idx, k_next) in enumerate(K_GRID)
                i = k_next - (1-delta)*k
                PROFIT_NOW = OP_PROFIT - ADJ_COST(model, i, k) - IRR(model, i)*i
                DISC_FUTURE_EXPECTED_PROFIT = beta*(sum(V_old[knxt_idx, :].*P_z[z_idx, :]))
                value = PROFIT_NOW + DISC_FUTURE_EXPECTED_PROFIT

                if value > best_value
                    V_new[k_idx, z_idx] = value
                    INVESTMENT_POLICY[k_idx, z_idx] = i
                    FUTURE_CAPITAL_POLICY[k_idx, z_idx] = k_next
                    ADJ_COST_POLICY[k_idx, z_idx] = ADJ_COST(model, i, k)
                    IRR_COST_POLICY[k_idx, z_idx] = IRR(model, i)*i
                    OPERATING_PROFIT_POLICY[k_idx, z_idx] = OP_PROFIT
                    TOTAL_PROFIT_POLICY[k_idx, z_idx] = PROFIT_NOW
                    best_value = value

                end
            end
        end
    end
    return V_new, INVESTMENT_POLICY, FUTURE_CAPITAL_POLICY, ADJ_COST_POLICY, IRR_COST_POLICY, OPERATING_PROFIT_POLICY, TOTAL_PROFIT_POLICY
end
    
function Solve_Model(model::ProjectParams, V_old, min_error, max_iter)
    error = 1
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
end #moduł