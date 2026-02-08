module Consumption_Savings_Model

using QuantEcon, Interpolations, LinearAlgebra, Parameters, Printf, Roots

export Consumption_Savings_Model

@with_kw struct Consumption_Savings_Model

    #Utiity function
    u = γ == 1 ? (c -> log(c)) : (c -> (c^(1 - γ) - 1) / (1 - γ))
    u_prime = γ == 1 ? (c -> 1/c) : (c -> c^(-γ))
    u_prime_inv = γ == 1 ? (y -> 1/y) : (y -> y^(-1/γ))

    γ = [2, 10] #setting two possible options for the risk aversion Parameters
    β = 0.99 # discount factor
    R = [1.010, 1.008] #interest rate -> R[1] when γ[1], and R[2] when γ[2], there might be some other smarter conditioning for Parameters

    #AR process for income
    ρ_z = 0.90 #income persistance
    σ_ϵ = 0.2*sqrt(1 - ρ^2) 
    N_z = 3  #Number of points in markov process
    mc_z = rouwenhorst(N_z, ρ_z, σ_ϵ, 0.0)
    λ_z = stationary_distributions(mc_z)[1]
    P_z = mc_z.p
    z_vec = exp.(mc_z.state_values) / sum(exp.(mc_z.state_values) .* λ_z)
    z_min = min(mc_z) #calculating z_min for the borrowing limit

    #asset grid 
    a_min = -0.6 .* z_min ./ (R .- 1) # minimum of the assets depend on the interest rate so it  is vectorized
    a_max = 500.0
    N_a = 100
    θ = 3
    ω = range(0, 1, length = N_a) # two different grids
    a_vec = a_min .+ (a_max - a_min) .* ω.^θ 

    L = sum(z_vec .* λ_z) #expected income
end
# how to change γ wihout weird dimensonality

end #ending module

