using Random, Distributions, StatsPlots, Statistics, Plots, NLopt

Random.seed!(2024)

#parameters
ρ = 0.9
p = 0.8
σ_L = 0.1
σ_H = 0.3
T = 500
θ = [ρ , p]
S = 100

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
observed_burnt = observed[101:end]
histogram(observed_burnt)
plot(observed_burnt)

#getting rid of logarythm to visualize
exponential_stochastic = @. exp(observed_burnt)
plot(exponential_stochastic) # looks similar to stock market 



# helper function for later

function moments_storage(data)
    T = length
    m_1 = std(data) 
    m_2 = cor(data[1 : end - 1], observed_burnt[2:end]) 
    m_3 = kurtosis(data[2:400] .- data[1:399]) 

    moments=[m_1,m_2,m_3]

    return moments
end

moments_storage(observed_burnt)

m = moments_storage(observed_burnt)

# displaying moments
m[1] # 0.27289655729182555
m[2] # 0.8539455367994457
m[3] # 3.443869988098574



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



function smm_objective(θ::Vector, gradient::Vector, observed_data::Vector, σ_L::Float64, σ_H::Float64, S::Int)
    Random.seed!(hash(θ)) 
    T_sim = length(observed_data)
    
    m1_list = []
    m2_list = []
    m3_list = []
    
    #list of moments 
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

    #means of those moments
    
    observed_moments=moments_storage(observed_data)
    #use the helper function to get the 'empirical' moments
    Q=(observed_moments[1]-sim_mean_m1)^2 +(observed_moments[2]-sim_mean_m2)^2 + (observed_moments[3]-sim_mean_m3)^2

    return Q 
end   


smm_objective(θ, [] ,observed_burnt, σ_L, σ_H,100)

#test
#plot(sims[2][:,1])
#histogram(sims[2][:,50])

# I did optimization and visual analysis

opt = NLopt.Opt(:LN_COBYLA, 2)

## Define the objective function:
NLopt.min_objective!(opt, (θ,gradient)->smm_objective(θ, gradient, observed_burnt, σ_L,σ_H, S))

## Define the lower bounds for the two parameters:
opt.lower_bounds = [0.5, 0.5] 
## Define the upper bounds for the two parameters:
opt.upper_bounds = [0.99, 0.95]   
## Define the stopping criteria:
opt.maxeval      = 2000
opt.xtol_rel     = 1e-10     
## Perform optimization on the object defined and the initial guess:
min_f, θ_optim, ret = NLopt.optimize(opt, [0.85, 0.7])

# comparison of optimized parameters vs real parameters - pretty close i guess but is it enough?
println(θ_optim) # θ_optim = [ρ =  0.8460851455325016,p = 0.7672396142192058]
println(θ)       # θ       = [ρ =  0.9, p = 0.8]

#simulating model with estimated theta
simulated_series = simulate_model(θ_optim, T, σ_L, σ_H)

# setting some trash theta to check if there is a difference for optimized theta
trash = [0.2, 0.5]
trash_series = simulate_model(trash, T, σ_L, σ_H)

# plotting series for comparison. We plot series with trash theta so we can compare if there is a value added of optimization
plot(observed_burnt[1:200], label = "observed series")
plot!(simulated_series[1:200], label = "simulated series")

#We plot series with trash theta so we can compare if there is a value added of optimization
plot(observed_burnt[1:200], label = "observed series",normalize=:pdf)
plot!(trash_series[1:200], label = "trash series", color = "green",normalize=:pdf)

#histograms of the data
histogram(observed_burnt, label = "observed series", normalize=:pdf, alpha = 0.6)
histogram!(simulated_series, label = " simulated series", normalize=:pdf, alpha = 0.5)
#histogram!(trash_series, label = "trash series", normalize=:pdf, alpha = 0.1)   # we can compare trash series distribution optionally


#calculating Lags 

Δ_observed  = observed_burnt[2 : end] .- observed_burnt[1 : end - 1]
Δ_simulated = simulated_series[2 : end] .- simulated_series[1 : end - 1]
Δ_trash     = trash_series[2:end] .- trash_series[1 : end-1] #also creating trash lag for comparison


#histograms of Δlog(yt)
histogram(Δ_observed, 
          label = "Observed series", 
          normalize = :pdf, 
          alpha = 0.6, 
          bins = 50, 
          title = "(Kurtosis)",
          xlabel = " (Δ log y)",
          ylabel = "(PDF)",
          legend = :topright)
histogram!(Δ_simulated, 
           label = "simulated series ρ^, p^", 
           normalize = :pdf, 
           alpha = 0.5,
           bins = 50)
histogram!(Δ_trash, 
           label = "Trash Series ( ρ=0.2, p=0.5)", 
           normalize = :pdf, 
           alpha = 0.2, 
           bins = 50, 
           color = "green")

#estimated parameters reproduce the key features of the observed data, because the shape is similar. 
#the behaviour of the process is similar