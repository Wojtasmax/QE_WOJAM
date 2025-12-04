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


#nie działą bo seed musi być w funkcji to wszystko wyzej powinno byc w funkcji 


# helper function for later to dział ale mozna to zorbic lepiej

function moments_storage(data)
    T = length
    m_1 = std(data) #0.2548018495221319
    m_2 = cor(data[1 : end - 1], observed_burnt[2:end]) #0.8334790574315395
    m_3 = kurtosis(data[2:400] .- data[1:399]) #3.3917448216937363

    moments=[m_1,m_2,m_3]

    return moments
end

moments_storage(observed_burnt)

m = moments_storage(observed_burnt)

# displaying moments
m[1]
m[2]
m[3]



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






#Olaf: zrobiłem ta funkcję 



function smm_objective(θ; observed_data, σ_L, σ_H, S)
    Random.seed!(hash(θ)) 
    T_sim = length(observed_data)
    
    m1_list = []
    m2_list = []
    m3_list = []
    
    #listy momentów dla wszystkich symulacji
    simulations = ones(T_sim,S)

    for i in 1:S
        simulations[:, i] =  simulate_model(θ, T_sim + 100, σ_L, σ_H)
        
        sim_m1 = std(simulations[:,i]) 
        sim_m2 = cor(simulations[1:end-1, i], simulations[2:end, i]) 
        sim_m3 = kurtosis(simulations[2:end,i] .- simulations[1:end-1,i])
        
        push!(m1_list, sim_m1)
        push!(m2_list, sim_m2)
        push!(m3_list,sim_m3)
    end

    sim_mean_m1 = mean(m1_list)
    sim_mean_m2 = mean(m2_list)
    sim_mean_m3 = mean(m3_list)

    #średnie z tych list momentów
    
    observed_moments=moments_storage(observed_data)
    #używam funkcji zeby wyciagnąć observed moments

    Q=(observed_moments[1]-sim_mean_m1)^2 +(observed_moments[2]-sim_mean_m2)^2 + (observed_moments[3]-sim_mean_m3)^2

    return (Q,simulations) #TODO later delete simulations from return, now to enable check
end   



sims=smm_objective(θ,observed_burnt,σ_L,σ_H,100)

#test
plot(sims[2][:,1])
histogram(sims[2][:,50])