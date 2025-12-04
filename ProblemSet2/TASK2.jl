using Random, Distributions, StatsPlots, Statistics, Plots, NLopt

Random.seed!(2024)

#parameters
ρ = 0.9
p = 0.8
σ_L = 0.1
σ_H = 0.3
T = 500

# creating epsilon function for parameters
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

#moments of observed_burnt

m_1 = std(observed_burnt) #0.2548018495221319
m_2 = cor(observed_burnt[1:399], observed_burnt[2:400]) #0.8334790574315395
m_3 = kurtosis(observed_burnt[2:400] .- observed_burnt[1:399]) #3.3917448216937363