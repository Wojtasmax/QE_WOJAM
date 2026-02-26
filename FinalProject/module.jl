module Engine

using QuantEcon, Statistics, Parameters, Interpolations, Optim

export ProjectParams

@with_kw struct ProjectParams
    #==== Parametry globalne ====#
    α::Float64 = 0.3
    v::Float64 = 0.6
    r::Float64 = 0.04
    δ::Float64 = 0.08
    w::Float64 = 1.0
    β::Float64 = 1.0 / (1.0 + r)
    
    #==== Parametry procesu Markowa ====#
    ρ::Float64 = 0.9
    σ_ε::Float64 = 0.12
    N_z::Int = 7  #to jest z czapy

    #==== Wyniki procesu Markowa ====#
    z_vec::Vector{Float64}
    P_z::Matrix{Float64}
    λ_z::Vector{Float64}
end


function ProjectParams(; α=0.3, v=0.6, r=0.04, δ=0.08, w=1.0, ρ=0.9, σ_ε=0.12, N_z=7)
    
    β = 1.0 / (1.0 + r)
    
   
    z_tilde = exp(-σ_ε^2 / (2 * (1 - ρ^2)))
    μ_logz = log(z_tilde)
    
    
    mc_z = rouwenhorst(N_z, ρ, σ_ε, μ_logz)
    P_z = mc_z.p
    λ_z = stationary_distributions(mc_z)[1]
    
    
    z_raw = exp.(mc_z.state_values)
    z_vec = z_raw ./ sum(z_raw .* λ_z)
    
    return ProjectParams(
        α=α, v=v, r=r, δ=δ, w=w, β=β, 
        ρ=ρ, σ_ε=σ_ε, N_z=N_z, 
        z_vec=z_vec, P_z=P_z, λ_z=λ_z
    )
end

end