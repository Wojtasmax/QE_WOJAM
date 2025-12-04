using Random, Distributions, StatsPlots, Statistics, Plots, NLopt

Random.seed!(2024)

#parameters
ρ = 0.9
p = 0.8
σ_L = 0.1
σ_H = 0.3
T = 500
θ = [ρ , p]

# creating epsilon function for parameters - gets the random epsilon
function ϵ(p = 0.8 , σ_L = 0.1, σ_H = 0.3)
    r = rand(Uniform(0,1))
    d1 = Normal(0, σ_L)
    d2 = Normal(0, σ_H) 
    if r < p 
        return rand(d1)
    else
        return rand(d2)
    end
end

# the way to call this function
ϵ()

#creating  observed vector
observed = Vector{Float64}(undef, 500)
#initial conditiion
observed[1] = 0
#generating time series
for t in 1 : (T - 1)
    observed[t + 1] = ρ * observed[t] + ϵ()
end

#visual analysis
histogram(observed)
plot(observed)

#burning the initial 100 observations
observed_burnt = observed[101:500]
histogram(observed_burnt)
plot(observed_burnt)

#getting rid of logarythm to visualize
exponential_stochastic = @. exp(observed_burnt)
plot(exponential_stochastic) # looks similar to stock market 

#moments of observed_burnt - with lag definition TODO seed is not working properly

m_1 = std(observed_burnt) #0.2548018495221319
m_2 = cor(observed_burnt[1:399], observed_burnt[2:400]) #0.8334790574315395
m_3 = kurtosis(observed_burnt[2:400] .- observed_burnt[1:399]) #3.3917448216937363


function simulate_model(θ, T , σ_L, σ_H)
    if T <= 100
        print("T must be higher than 100")
    end

    simulated = Vector{Float64}(undef, T)
    simulated[1] = 0

    for t in 1 : (T - 1)
    simulated[t + 1] = θ[1] * simulated[t] + ϵ(θ[2], σ_L, σ_H)
    end
    return simulated[101 : T]
end

# safety check 
test = simulate_model(θ,500, σ_L, σ_H)
plot(test)
histogram(test)

std(test)

#TODO create seed for this function
function smm_objective(θ, observed_data, σ_L, σ_H, S)
    T = length(observed_data)
    
    #moments from data
    m_1 = std(observed_data)
    m_2 = cor(observed_data[1: length(observed_data - 1)], observed_data[2: length(observed_data)])
    m_3 = kurtosis(observed_data[2: length(observed_data)] .- observed_data[1: length(observed_data - 1)])

    simulations = Vector{Float64}(undef, S) # we need to store simulations in the matrix or vector of vectors?
    for i in 1:S 
        simulations[i] =  simulate_model(θ, T, σ_L, σ_H)
    end
    
    Simulation_Matrix =  
end    



