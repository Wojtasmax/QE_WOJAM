using Distributions, Plots, NLsolve, Roots


function model_parameters(γ;
    A=1.0, # TFP
    δ=0.1, # depreciation rate
    β=0.96, # discount factor
    α=0.33 # capital share
    )
    khat = (1/(A*α) * (1/β - 1 + δ))^(1/(α-1)) # steady state capital
    chat = A*khat^α - δ*khat # steady state consumption
    k0 = 0.5*khat # initial capital
    return (; γ, A, δ, β, α, khat, chat, k0)
end

function transition_equation(res, x, T, parameters)

    A  = parameters.A
    δ  = parameters.δ
    β  = parameters.β
    α  = parameters.α
    γ  = parameters.γ
    k0 = parameters.k0
    chat = parameters.chat

    res = zeros(2*T+1) # vector for residuals

    for t in 1:T
        res[t] = x[t]^(-γ) - (β*x[t+1]^(-γ) * (α*A*x[T+t+1]^(α-1) + 1 - δ)) # difference between consumption (in fact - the inverse of utility) from Euler equation and our data
    end

    for t in 2:T
        res[T+t+1] = x[T+t+1] - ((1 - δ)*x[T+t] + A*x[T+t]^α - x[t-1]) 
    end

    res[T+2] = x[T+2]-((1 - δ)*k0 + A*k0^α - x[1])
    res[T+1] = x[T+1]-chat
    return res
end

function simulate_path(c0, T, parameters)
    A = parameters.A
    δ = parameters.δ
    β = parameters.β
    α = parameters.α
    γ = parameters.γ

    c = zeros(T+1)
    k = zeros(T+1)

    c[1] = c0
    k[1] = parameters.k0

    for t in 1:T

        # capital evolution
        k[t+1] = (1 - δ)*k[t] + A*k[t]^α - c[t]

        if k[t+1] <= 0
            # return huge numbers so residual is very wrong → solver avoids this path
            return fill(10, T+1), fill(10, T+1)
        end

        # Euler equation
        c[t+1] = (β * c[t]^γ * (α*A*k[t+1]^(α-1) + 1 - δ))^(1/γ)

    end

    return c, k
end

function shooting_residual(z, T, parameters)
    c0 = exp(z)  
    cT, kT = simulate_path(c0, T, parameters)
    return cT[end] - parameters.chat
end


function find_initial_consumption(T, parameters)
    f(z) = shooting_residual(z, T, parameters)

    # initial bracket (in log space)
    z_low = log(0.1 * parameters.chat)
    z_high = log(2.0 * parameters.chat)

    # Expand bracket until signs differ (simple heuristic)
    for i in 1:20
        if f(z_low) * f(z_high) < 0
            break
        end
        z_low -= 0.5
        z_high += 0.5
    end

    # If bracket found, use Brent; otherwise fall back to find_zero with Order1 (Newton)
    if f(z_low) * f(z_high) < 0
        z_star = find_zero(f, (z_low, z_high), Bisection())
    else
        # fallback: Newton with initial guess log(chat)
        z_star = find_zero(f, log(parameters.chat), Order1())
    end

    return exp(z_star)
end

my_model = model_parameters(0.5);
find_initial_consumption(50, my_model);

c_1, k_1 = simulate_path(find_initial_consumption(50, my_model),50,my_model);

