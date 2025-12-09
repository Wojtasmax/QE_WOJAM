using Distributions, Plots


function steady_state(;
    A=1.0, # TFP
    δ=0.1, #depreciation rate
    β=0.96, # discount factor
    α=0.33 # capital share
    )
    khat = (1/(A*α) * (1/β - 1 + δ))^(1/(α-1))
    chat = A*khat^α - δ*khat
    k_begin = 0.5*khat
    return (khat, chat, k_begin)
end

function transition_equation()
end
