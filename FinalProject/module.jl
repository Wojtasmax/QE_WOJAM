module Engine

using QuantEcon, Statistics, Parameters, Interpolations, Optim

export ProjectParams

@with_kw struct ProjectParams
    #==== Parametry globalne ====#
    α= 0.3
    v= 0.6
    r= 0.04
    δ= 0.08
    w= 1.0
    β = 1.0 / (1.0 + r)
    
    #==== Parametry procesu Markowa ====#
    ρ::Float64 = 0.9
    σ_ε::Float64 = 0.12
    N_z::Int = 7  #to jest z czapy

    #==== Wyniki procesu Markowa ====#
    z_vec::Vector{Float64}=Float64[]
    P_z::Matrix{Float64}=zeros(1,1)
    λ_z::Vector{Float64}=Float64[]


    #===Adjustment Cost====#
    #losowe wartosci    
    γ = 0.05          
    F = 0.01          #fixed cost
    ps= 0.80         
end


function ProjectParams(; α=0.3, v=0.6, r=0.04, δ=0.08, w=1.0, ρ=0.9, σ_ε=0.12, N_z=7)
    @unpack α,v,r,δ,w,ρ,σ_ε,N_z

    β = 1.0 / (1.0 + r)
    
   
    z̃ = exp(-σ_ε^2 / (2 * (1 - ρ^2)))
    μ_logz = log(z̃)
    
    
    mc_z = rouwenhorst(N_z, ρ, σ_ε, μ_logz)
    P_z = mc_z.p
    λ_z = stationary_distributions(mc_z)[1]
    
    
    z_raw = exp.(mc_z.state_values)
    z_vec = z_raw ./ sum(z_raw .* λ_z)
    
    return ProjectParams(
        α=α, v=v, r=r, δ=δ, w=w, β=β, 
        ρ=ρ, σ_ε=σ_ε, N_z=N_z, 
        z_vec=z_vec, P_z=P_z, λ_z=λ_z,
        γ=γ, F=F, ps=ps
    )
end
function ADJ_COST(model::ProjectParams,i,k, )
    @unpack F,γ =model
    if i !=0
        return ((γ/2)*(i/k)^2)*k+F*k
    else
        return 0.0
    end
end


function IRR(model::ProjectParams,i)
    @unpack ps = model
    if i >=0
        return 1.0
    else
        return ps
    end
end


function OPTIMAL_H(model::ProjectParams,k,z) #jak cos sie tu bedzie jebac to zmienić tak zeby ulamke w potedze byl dodatni
    @unpack w,v,α = model
    return (w/(v*z*k^α))^(1/(v-1))
end


function OPERATING_PROFIT(model::ProjectParams,k,z)  
    @unpack α,v,w = model
    h=OPTIMAL_H(model::ProjectParams,k)
    return z*k^α*h^v-w*h
end
end #moduł