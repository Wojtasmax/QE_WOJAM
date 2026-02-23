using Plots, Printf, Parameters, Statistics, Random, Interpolations, LinearAlgebra
include("module.jl")

using .Consumption_Savings_Model: Consumption_Savings, vfi_interp_transformed, get_consumption_matrix

# Set seed for reproducibility
Random.seed!(1234)

# ============================================================================
# SETUP MODEL
# ============================================================================

model = Consumption_Savings(γ = 2.0, R = 1.010)

# Solve model using TRANSFORMED VFI
println("\nSolving model with transformed VFI...")
@time v, σ, iter, err = vfi_interp_transformed(model)

println("\nConvergence results:")
println("  Iterations: $iter")
println("  Error:      $err")

# Get consumption matrix
consumption = get_consumption_matrix(model, σ)

# ============================================================================
# MPC CALCULATION FUNCTION
# ============================================================================

function calculate_mpc(model, σ, δ)
    @unpack N_a, N_z, a_vec, z_vec, R = model
    
    mpc = zeros(N_a, N_z)
    
    # Create interpolators for policy function
    σ_interps = [LinearInterpolation(a_vec, σ[:, iz], extrapolation_bc=Line()) for iz in 1:N_z]
    
    for iz in 1:N_z
        z = z_vec[iz]
        σ_interp = σ_interps[iz]
        
        for ia in 1:N_a
            a = a_vec[ia]
            
            a_next = σ_interp(a)
            c_current = R * a + z - a_next
            
            a_plus_delta = a + δ
            a_next_plus = σ_interp(a_plus_delta)
            c_plus_delta = R * a_plus_delta + z - a_next_plus
            
            if c_current > 0 && c_plus_delta > 0
                mpc[ia, iz] = (c_plus_delta - c_current) / δ
            else
                mpc[ia, iz] = NaN
            end
        end
    end
    
    return mpc
end

# ============================================================================
# CALCULATE MPC FOR MULTIPLE δ VALUES
# ============================================================================

δ_values = [0.01, 0.1, 0.5, 1.0, 2.0]

# Store MPC results for each δ
mpc_dict = Dict{Float64, Matrix{Float64}}()

for δ in δ_values

    mpc_dict[δ] = calculate_mpc(model, σ, δ)
end

# ============================================================================
# PLOTTING MPC FOR EACH δ VALUE
# ============================================================================

# Income state indices
iz_mid = model.N_z ÷ 2 + 1
iz_low = 1
iz_high = model.N_z

# Plot MPC(a) for 3 income levels, separately for each δ
for δ in δ_values
    p = plot(model.a_vec, mpc_dict[δ][:, iz_low],
             label="z=$(round(model.z_vec[iz_low], digits=3))",
             linewidth=2, xlabel="Assets (a)", ylabel="MPC",
             title="MPC (δ=$δ, γ=2.0, R=1.010)",
             legend=:topright, ylims=(0, 1.2))
    
    plot!(p, model.a_vec, mpc_dict[δ][:, iz_mid],
          label="z=$(round(model.z_vec[iz_mid], digits=3))",
          linewidth=2, linestyle=:dash)
    
    plot!(p, model.a_vec, mpc_dict[δ][:, iz_high],
          label="z=$(round(model.z_vec[iz_high], digits=3))",
          linewidth=2, linestyle=:dot)
    
    # Add reference line at MPC = 1
    plot!(p, [model.a_min, model.a_max], [1.0, 1.0],
          label="MPC = 1", linewidth=1, color=:black, linestyle=:dash, alpha=0.5)
    
    display(p)
end

# ============================================================================
# SOLVING FOR STATIONARY DISTRIBUTION AND AGGREGATE CONSUMPTION
# ============================================================================

function compute_stationary_distribution(model, σ; maxiter=20000, tol=1e-7)
    @unpack N_a, N_z, a_vec, z_vec, P_z, λ_z = model
    
    # Initialize distribution: start with ergodic distribution of z, uniform over a
    λ = ones(N_a, N_z) .* λ_z' / N_a  # Normalize so sum = 1
    
    # Create interpolators for policy function
    σ_interps = [LinearInterpolation(a_vec, σ[:, iz], extrapolation_bc=Line()) for iz in 1:N_z]
    
    for iter in 1:maxiter
        λ_new = zeros(N_a, N_z)
        
        # For each current state (a, z)
        for iz in 1:N_z
            σ_interp = σ_interps[iz]
            
            for ia in 1:N_a
                a = a_vec[ia]
                a_next = σ_interp(a)  # Where this agent moves to
                
                # Find which grid points a_next falls between
                # and distribute mass proportionally (linear interpolation)
                if a_next <= a_vec[1]
                    ia_next = 1
                    weight_lower = 1.0
                    weight_upper = 0.0
                elseif a_next >= a_vec[end]
                    ia_next = N_a - 1
                    weight_lower = 0.0
                    weight_upper = 1.0
                else
                    # Find bracketing indices
                    ia_next = searchsortedlast(a_vec, a_next)
                    ia_next = min(ia_next, N_a - 1)  # Ensure we don't go out of bounds
                    
                    # Linear interpolation weights
                    a_lower = a_vec[ia_next]
                    a_upper = a_vec[ia_next + 1]
                    weight_upper = (a_next - a_lower) / (a_upper - a_lower)
                    weight_lower = 1.0 - weight_upper
                end
                
                # Distribute mass to next period
                for iz_next in 1:N_z
                    prob = P_z[iz, iz_next]  # Probability of transitioning to iz_next
                    
                    if ia_next <= N_a - 1
                        λ_new[ia_next, iz_next] += weight_lower * prob * λ[ia, iz]
                        λ_new[ia_next + 1, iz_next] += weight_upper * prob * λ[ia, iz]
                    else
                        λ_new[ia_next, iz_next] += prob * λ[ia, iz]
                    end
                end
            end
        end
        
        # Check convergence
        diff = maximum(abs.(λ_new - λ))
        
        if diff < tol
            println("  Converged in $iter iterations, max difference: $diff")
            return λ_new
        end
        
        λ = λ_new
    end
    
    println("  Warning: Did not converge after $maxiter iterations")
    return λ
end

# Compute stationary distributions of assets and consumption
λ_star = compute_stationary_distribution(model, σ);

function compute_aggregate_consumption(consumption, λ_star)
    return sum(consumption .* λ_star)
end

C_star = compute_aggregate_consumption(consumption, λ_star);

#Transfer size
Δ = 0.05 * C_star;

# ============================================================================
# GIVING TRANSFERS AND SHIFTING DISTRIBUTIONS
# ============================================================================

function shift_distribution(model, λ_star, Δ)
    @unpack N_a, N_z, a_vec = model
    
    λ_0 = zeros(N_a, N_z)
    
    for iz in 1:N_z
        for ia in 1:N_a
            a = a_vec[ia]
            a_shifted = a + Δ  # New asset position after transfer
            
            # Find where a_shifted falls on the grid
            if a_shifted <= a_vec[1]
                # Below grid, put all mass on first point
                λ_0[1, iz] += λ_star[ia, iz]
            elseif a_shifted >= a_vec[end]
                # Above grid, put all mass on last point
                λ_0[N_a, iz] += λ_star[ia, iz]
            else
                # Interpolate between grid points
                ia_next = searchsortedlast(a_vec, a_shifted)
                ia_next = min(ia_next, N_a - 1)
                
                a_lower = a_vec[ia_next]
                a_upper = a_vec[ia_next + 1]
                weight_upper = (a_shifted - a_lower) / (a_upper - a_lower)
                weight_lower = 1.0 - weight_upper
                
                # Distribute mass
                λ_0[ia_next, iz] += weight_lower * λ_star[ia, iz]
                λ_0[ia_next + 1, iz] += weight_upper * λ_star[ia, iz]
            end
        end
    end
    
    return λ_0
end

#Shifting distribution
λ_0 = shift_distribution(model, λ_star, Δ);

#Solving for new consumption
C_0 = compute_aggregate_consumption(consumption, λ_0);

#Calculating aggregate MPC
aggregate_MPC = (C_0 - C_star) / Δ;

# ============================================================================
# IMPULSE RESPONSE ANALYSIS
# ============================================================================

function forward_iterate(model, σ, λ_t)
    @unpack N_a, N_z, a_vec, P_z = model
    
    λ_next = zeros(N_a, N_z)
    
    # Create interpolators for policy function
    σ_interps = [LinearInterpolation(a_vec, σ[:, iz], extrapolation_bc=Line()) for iz in 1:N_z]
    
    # For each current state (a, z)
    for iz in 1:N_z
        σ_interp = σ_interps[iz]
        
        for ia in 1:N_a
            a = a_vec[ia]
            a_next = σ_interp(a)  # Where this agent moves to
            
            # Find which grid points a_next falls between
            if a_next <= a_vec[1]
                ia_next = 1
                weight_lower = 1.0
                weight_upper = 0.0
            elseif a_next >= a_vec[end]
                ia_next = N_a - 1
                weight_lower = 0.0
                weight_upper = 1.0
            else
                # Find bracketing indices
                ia_next = searchsortedlast(a_vec, a_next)
                ia_next = min(ia_next, N_a - 1)
                
                # Linear interpolation weights
                a_lower = a_vec[ia_next]
                a_upper = a_vec[ia_next + 1]
                weight_upper = (a_next - a_lower) / (a_upper - a_lower)
                weight_lower = 1.0 - weight_upper
            end
            
            # Distribute mass to next period
            for iz_next in 1:N_z
                prob = P_z[iz, iz_next]
                
                if ia_next <= N_a - 1
                    λ_next[ia_next, iz_next] += weight_lower * prob * λ_t[ia, iz]
                    λ_next[ia_next + 1, iz_next] += weight_upper * prob * λ_t[ia, iz]
                else
                    λ_next[ia_next, iz_next] += prob * λ_t[ia, iz]
                end
            end
        end
    end
    
    return λ_next
end

# ============================================================================
# FORWARD SIMULATION: t = 0, 1, ..., 50
# ============================================================================

T_max = 50
λ_path = Vector{Matrix{Float64}}(undef, T_max + 1)
C_path = Vector{Float64}(undef, T_max + 1)

# Initial condition: post-transfer distribution
λ_path[1] = λ_0  # t = 0
C_path[1] = C_0

# Forward iteration
for t in 1:T_max
    λ_path[t + 1] = forward_iterate(model, σ, λ_path[t])
    C_path[t + 1] = compute_aggregate_consumption(consumption, λ_path[t + 1])
end

# ============================================================================
# 1. PLOT IMPULSE RESPONSE
# ============================================================================

# Compute impulse response: (C^t - C*) / C*
time_grid = 0:T_max
impulse_response = [(C_path[t+1] - C_star) / C_star for t in 0:T_max]

p_impulse = plot(time_grid, impulse_response,
                 label="Impulse Response",
                 linewidth=2, marker=:circle, markersize=3,
                 xlabel="Time (t)", ylabel="(Cᵗ - C*) / C*",
                 title="Consumption Impulse Response to Transfer",
                 legend=:topright, grid=true)
hline!(p_impulse, [0.0], label="Steady State", linestyle=:dash, color=:black, alpha=0.5)
display(p_impulse)

# ============================================================================
# 2. COMPUTE CUMULATIVE MPCs
# ============================================================================

function cumulative_mpc(C_path, C_star, Δ, H)
    if H == 0
        return 0.0
    end
    sum_consumption_increase = sum(C_path[t+1] - C_star for t in 0:H-1)
    return sum_consumption_increase / Δ
end

# Horizons to compute
H_values = [0, 1, 4, 8, 12, 20]

for H in H_values
    mpc_H = cumulative_mpc(C_path, C_star, Δ, H)
end

# ============================================================================
# 3. FRACTION OF TRANSFER SPENT BY H = 12
# ============================================================================

# Fraction spent = MPC^H (since transfer is Δ and we compute cumulative consumption increase)
H_report = 12
mpc_12 = cumulative_mpc(C_path, C_star, Δ, H_report)
fraction_spent = mpc_12

# ============================================================================
# PLOT CUMULATIVE MPC OVER TIME
# ============================================================================

# Compute cumulative MPC for all horizons
cumulative_mpc_path = [cumulative_mpc(C_path, C_star, Δ, H) for H in 1:T_max+1]

p_cumulative = plot(1:T_max+1, cumulative_mpc_path,
                    label="Cumulative MPC",
                    linewidth=2, marker=:circle, markersize=2,
                    xlabel="Horizon H", ylabel="MPC^H",
                    title="Cumulative MPC over Horizon",
                    legend=:bottomright, grid=true)

# Add markers for specific horizons
for H in [1, 4, 8, 12, 20]
    if H <= T_max + 1
        scatter!(p_cumulative, [H], [cumulative_mpc_path[H]], 
                marker=:star, markersize=8, color=:red, label="")
    end
end

hline!(p_cumulative, [1.0], label="MPC = 1 (full spending)", linestyle=:dash, color=:black, alpha=0.5)
display(p_cumulative)

# ============================================================================
# ADDITIONAL ANALYSIS: CONSUMPTION INCREASE OVER TIME
# ============================================================================

# Absolute consumption increase
consumption_increase = [C_path[t+1] - C_star for t in 0:T_max]

p_cons_increase = plot(time_grid, consumption_increase,
                       label="Cᵗ - C*",
                       linewidth=2, marker=:circle, markersize=3,
                       xlabel="Time (t)", ylabel="Cᵗ - C*",
                       title="Aggregate Consumption Increase over Time",
                       legend=:topright, grid=true)
hline!(p_cons_increase, [0.0], label="Steady State", linestyle=:dash, color=:black, alpha=0.5)
display(p_cons_increase)

# ============================================================================
# SUMMARY TABLE
# ============================================================================

println("\n" * "="^70)
println("SUMMARY OF RESULTS")
println("="^70)
println("\nSteady State:")
@printf("  C* = %.6f\n", C_star)
println("\nTransfer:")
@printf("  Δ = 0.05 × C* = %.6f\n", Δ)
println("\nImmediate Response (t=0):")
@printf("  C⁰ = %.6f\n", C_0)
@printf("  Impact MPC = (C⁰ - C*) / Δ = %.4f\n", aggregate_MPC)
println("\nCumulative MPCs:")
for H in H_values
    if H <= T_max
        mpc_H = cumulative_mpc(C_path, C_star, Δ, H)
        @printf("  MPC^%d  = %.4f (%.2f%% of transfer)\n", H, mpc_H, mpc_H * 100)
    end
end
println("\nFraction Spent by H=12:")
@printf("  %.2f%% of transfer\n", fraction_spent * 100)