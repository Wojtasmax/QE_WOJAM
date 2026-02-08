using Distributions, LinearAlgebra,Plots, Random, QuantEcon

Random.seed!(1234)

P = [0.6 0.3 0.1; 0.2 0.6 0.2; 0.1 0.3 0.6]

periods = 1000


## Define a function that simulates a Markov chain
function mc_sample_path(P; init_x = 1, sample_size = 100)
    @assert size(P)[1] == size(P)[2] # square required
    N = size(P)[1] # number of states

    # Translate rows of transition matrix P
    # into a vector of distributions of discrete RV 
    dists = [Categorical(P[i, :]) for i in 1:N]
    
    # Setup the simulation
    X = Vector{Int64}(undef, sample_size) # allocate memory
    X[1] = init_x # set the initial state

    for t in 2:sample_size
        previous_state_value = X[t-1]
        P_Xt = dists[previous_state_value] # appropriate distribution  
        X[t] = rand(P_Xt) # draw new value
    end
    return X
end


sample_path_initSR = mc_sample_path(P, init_x = 3, sample_size = periods);
Z_t = [1.0 , 2.0 , 3.0]; 
X_t = [0.0 1.0 2.0 3.0 4.0 5.0];


time_series = Z_t[sample_path_initSR];

plot(time_series,xlabel = "time",ylabel = "Z_t", label=false)

time = [i for i in 0:5]

time_series[1]


function σ(X, Z)
    if Z == 1
        return 0
    end
    if Z == 2
        return X
    end
    if Z == 3 && X <= 4 
        return min(X + 1, 5)
    end
    if Z == 3 && X == 5
        return 3
    end
    return X 
end


state_space = Matrix{Int64}(undef, 18, 3)

X_vals = 0:5
Z_vals = 1:3
n_X = length(X_vals)
n_Z = length(Z_vals)
n_states = n_X * n_Z


state_list = [(x, z) for x in X_vals for z in Z_vals]

Q = zeros(n_states, n_states)

for i in 1:n_states
    x, z = state_list[i]
    
    x_next = σ(x, z)
    
    for j in 1:n_states
        x_prime, z_prime = state_list[j]
        
        if x_prime == x_next
            Q[i, j] = P[z, z_prime]
        else
            Q[i, j] = 0.0
        end
    end
end

println("Size of Q: ", size(Q))


mc = MarkovChain(Q)
ψ_star = stationary_distributions(mc)[1]


ψ_matrix = reshape(ψ_star, (n_Z, n_X))


ψ_X = vec(sum(ψ_matrix, dims=1)) 


ψ_Z = vec(sum(ψ_matrix, dims=2))


E_X = sum(X_vals .* ψ_X)
println("Mean of X: ", round(E_X, digits=4))

E_X_given_Z = zeros(n_Z)

for z_idx in 1:n_Z
    joint_probs_given_z = ψ_matrix[z_idx, :]
    

    cond_probs = joint_probs_given_z / ψ_Z[z_idx]
    

    E_X_given_Z[z_idx] = sum(X_vals .* cond_probs)
end

println("Conditional Means E[X | Z]: ", round.(E_X_given_Z, digits=4))


p1 = bar(X_vals, ψ_X, 
    title="Marginal Distribution of X", 
    xlabel="X", ylabel="Probability", 
    legend=false, color=:blue
)

p2 = bar(Z_vals, ψ_Z, 
    title="Marginal Distribution of Z", 
    xlabel="Z", ylabel="Probability", 
    legend=false, color=:red
)

p3 = bar(Z_vals, E_X_given_Z, 
    title="Conditional Mean E[X | Z]", 
    xlabel="Z", ylabel="Expected X", 
    legend=false, color=:green
)

final_plot = plot(p1, p2, p3, layout=(3,1), size=(600, 800))

display(final_plot)