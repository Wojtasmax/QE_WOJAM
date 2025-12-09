using Distributions, Plots, NLsolve


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

function solve_transition_path(T, parameters)

    chat = parameters.chat
    khat = parameters.khat
    k0 = parameters.k0

    c_guess = fill(chat, T+1)
    k_guess = range(k0, khat, length=T+1)[2:end]
    x0 = vcat(c_guess, k_guess)

    sol = nlsolve((res, x) -> transition_equation(res, x, T, parameters), x0,  method=:newton)

    if sol.converged
        println("Converged")
    else 
        println("Not converged")
    end

    return sol
end


my_model = model_parameters(0.5)
solve_transition_path(200, my_model)