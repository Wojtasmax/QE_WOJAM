
using Plots, Printf, Parameters
include("module.jl")

using .Consumption_Savings_Model: Consumption_Savings, T_interp, vfi_interp, create_initial_guess, get_consumption_matrix, euler_residuals

model_low_aversion  = Consumption_Savings(γ = 2.0 , R = 1.010)
model_high_aversion = Consumption_Savings(γ = 10.0 , R = 1.008)

println("Creating initial guess for all aversion cases...")
v_init = create_initial_guess(model_low_aversion)


println("Solving low risk aversion case...")
@time v_low, σ_low, iter_low, err_low = vfi_interp(model_low_aversion, v_init=v_init)
#VFI with interpolation converged in 1000 iterations, error: 0.0003244734177545071
#3.476393 seconds (13.96 M allocations: 458.864 MiB, 9.36% gc time, 35.51% compilation time)

println("Solving high risk aversion case...")
@time v_high, σ_high, iter_high, err_high = vfi_interp(model_high_aversion, v_init=v_init)
#VFI with interpolation converged in 1000 iterations, error: 0.01000022955375717
#2.450367 seconds (11.08 M allocations: 308.354 MiB, 2.41% gc time)


# Indices for different income states
iz_mid = model_low_aversion.N_z ÷ 2 + 1
iz_low = 1
iz_high = model_low_aversion.N_z

p1 = plot(model_low_aversion.a_vec, v_low[:, iz_mid], label="Median income", linewidth=2,
          xlabel="a", ylabel="v(a,z)", title="Value Function (z=$(round(model_low_aversion.z_vec[iz_mid], digits=3)))")
plot!(p1, model_low_aversion.a_vec, v_low[:, iz_low], label="low income", linewidth=2, linestyle=:dash)
plot!(p1, model_low_aversion.a_vec, v_low[:, iz_high], label="High income", linewidth=2, linestyle=:dot)
display(p1)
#low aversion case: value function is increasing in assets, and higher for higher income states. The gap between income states is larger at low asset levels, and narrows as assets increase.
p1b= plot(model_low_aversion.a_vec, v_low[:, iz_low], label="z=$(round(model_low_aversion.z_vec[iz_low], digits=3))", 
           linewidth=2, xlabel="a", ylabel="v(a,z)", title="Value Function (Different z)")
plot!(p1b, model_low_aversion.a_vec, v_low[:, iz_mid], label="z=$(round(model_low_aversion.z_vec[iz_mid], digits=3))", 
      linewidth=2, linestyle=:dash)
plot!(p1b, model_low_aversion.a_vec, v_low[:, iz_high], label="z=$(round(model_low_aversion.z_vec[iz_high], digits=3))", 
      linewidth=2, linestyle=:dot)
plot(p1b)

#high aversion case: value function is much more concave, and the gap between income states is much larger, especially at low asset levels. The value of having more assets is much higher for low income states, as it provides a buffer against bad shocks.
p1c = plot(model_high_aversion.a_vec, v_high[:, iz_low], label="z=$(round(model_high_aversion.z_vec[iz_low], digits=3))", 
           linewidth=2, xlabel="a", ylabel="v(a,z)", title="Value Function (Different z)")
plot!(p1c, model_high_aversion.a_vec, v_high[:, iz_mid], label="z=$(round(model_high_aversion.z_vec[iz_mid], digits=3))", 
      linewidth=2, linestyle=:dash)
plot!(p1c, model_high_aversion.a_vec, v_high[:, iz_high], label="z=$(round(model_high_aversion.z_vec[iz_high], digits=3))", 
      linewidth=2, linestyle=:dot)

# Zoom in near budget constraint (a_min) — use a small fraction of the grid for tighter zoom
n_zoom_frac = 0.02  # fraction of the asset grid to show (2%)
idx_low = max(2, ceil(Int, n_zoom_frac * length(model_low_aversion.a_vec)))
zoom_range_low = model_low_aversion.a_vec[idx_low]
p1b_zoom = plot(model_low_aversion.a_vec, v_low[:, iz_low], label="z=$(round(model_low_aversion.z_vec[iz_low], digits=3))", 
                 linewidth=2, xlabel="a", ylabel="v(a,z)", title="Value Function near a_min (Low Aversion)")
plot!(p1b_zoom, model_low_aversion.a_vec, v_low[:, iz_mid], label="z=$(round(model_low_aversion.z_vec[iz_mid], digits=3))", 
      linewidth=2, linestyle=:dash)
plot!(p1b_zoom, model_low_aversion.a_vec, v_low[:, iz_high], label="z=$(round(model_low_aversion.z_vec[iz_high], digits=3))", 
      linewidth=2, linestyle=:dot)
xlims!(p1b_zoom, model_low_aversion.a_min, zoom_range_low)

idx_high = max(2, ceil(Int, n_zoom_frac * length(model_high_aversion.a_vec)))
zoom_range_high = model_high_aversion.a_vec[idx_high]
p1c_zoom = plot(model_high_aversion.a_vec, v_high[:, iz_low], label="z=$(round(model_high_aversion.z_vec[iz_low], digits=3))", 
                 linewidth=2, xlabel="a", ylabel="v(a,z)", title="Value Function near a_min (High Aversion)")
plot!(p1c_zoom, model_high_aversion.a_vec, v_high[:, iz_mid], label="z=$(round(model_high_aversion.z_vec[iz_mid], digits=3))", 
      linewidth=2, linestyle=:dash)
plot!(p1c_zoom, model_high_aversion.a_vec, v_high[:, iz_high], label="z=$(round(model_high_aversion.z_vec[iz_high], digits=3))", 
      linewidth=2, linestyle=:dot)
xlims!(p1c_zoom, model_high_aversion.a_min, zoom_range_high)

# Medium zoom (10% of grid)
n_medium_frac = 0.10  # fraction of the asset grid to show (10%)
idx_low_med = max(2, ceil(Int, n_medium_frac * length(model_low_aversion.a_vec)))
zoom_med_low = model_low_aversion.a_vec[idx_low_med]
idx_high_med = max(2, ceil(Int, n_medium_frac * length(model_high_aversion.a_vec)))
zoom_med_high = model_high_aversion.a_vec[idx_high_med]

# Value function - medium zoom
p1b_med = plot(model_low_aversion.a_vec, v_low[:, iz_low], label="z=$(round(model_low_aversion.z_vec[iz_low], digits=3))",
                linewidth=2, xlabel="a", ylabel="v(a,z)", title="Value Function medium zoom (Low Aversion)")
plot!(p1b_med, model_low_aversion.a_vec, v_low[:, iz_mid], label="z=$(round(model_low_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p1b_med, model_low_aversion.a_vec, v_low[:, iz_high], label="z=$(round(model_low_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
xlims!(p1b_med, model_low_aversion.a_min, zoom_med_low)
display(p1b_med)

p1c_med = plot(model_high_aversion.a_vec, v_high[:, iz_low], label="z=$(round(model_high_aversion.z_vec[iz_low], digits=3))",
                linewidth=2, xlabel="a", ylabel="v(a,z)", title="Value Function medium zoom (High Aversion)")
plot!(p1c_med, model_high_aversion.a_vec, v_high[:, iz_mid], label="z=$(round(model_high_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p1c_med, model_high_aversion.a_vec, v_high[:, iz_high], label="z=$(round(model_high_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
xlims!(p1c_med, model_high_aversion.a_min, zoom_med_high)
display(p1c_med)

# Policy Functions - medium zoom
p2b_med = plot(model_low_aversion.a_vec, σ_low[:, iz_low], label="z=$(round(model_low_aversion.z_vec[iz_low], digits=3))",
                linewidth=2, xlabel="a", ylabel="σ(a,z)", title="Policy Function medium zoom (Low Aversion)")
plot!(p2b_med, model_low_aversion.a_vec, σ_low[:, iz_mid], label="z=$(round(model_low_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p2b_med, model_low_aversion.a_vec, σ_low[:, iz_high], label="z=$(round(model_low_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
xlims!(p2b_med, model_low_aversion.a_min, zoom_med_low)
display(p2b_med)

p2c_med = plot(model_high_aversion.a_vec, σ_high[:, iz_low], label="z=$(round(model_high_aversion.z_vec[iz_low], digits=3))",
                linewidth=2, xlabel="a", ylabel="σ(a,z)", title="Policy Function medium zoom (High Aversion)")
plot!(p2c_med, model_high_aversion.a_vec, σ_high[:, iz_mid], label="z=$(round(model_high_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p2c_med, model_high_aversion.a_vec, σ_high[:, iz_high], label="z=$(round(model_high_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
xlims!(p2c_med, model_high_aversion.a_min, zoom_med_high)
display(p2c_med)

# Consumption - medium zoom
p3b_med = plot(model_low_aversion.a_vec, consumption_low_aversion[:, iz_low], label="z=$(round(model_low_aversion.z_vec[iz_low], digits=3))",
                linewidth=2, xlabel="a", ylabel="c(a,z)", title="Consumption medium zoom (Low Aversion)")
plot!(p3b_med, model_low_aversion.a_vec, consumption_low_aversion[:, iz_mid], label="z=$(round(model_low_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p3b_med, model_low_aversion.a_vec, consumption_low_aversion[:, iz_high], label="z=$(round(model_low_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
xlims!(p3b_med, model_low_aversion.a_min, zoom_med_low)
display(p3b_med)

p3c_med = plot(model_high_aversion.a_vec, consumption_high_aversion[:, iz_low], label="z=$(round(model_high_aversion.z_vec[iz_low], digits=3))",
                linewidth=2, xlabel="a", ylabel="c(a,z)", title="Consumption medium zoom (High Aversion)")
plot!(p3c_med, model_high_aversion.a_vec, consumption_high_aversion[:, iz_mid], label="z=$(round(model_high_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p3c_med, model_high_aversion.a_vec, consumption_high_aversion[:, iz_high], label="z=$(round(model_high_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
xlims!(p3c_med, model_high_aversion.a_min, zoom_med_high)
display(p3c_med)

# Policy Functions - Full Range
p2b = plot(model_low_aversion.a_vec, σ_low[:, iz_low], label="z=$(round(model_low_aversion.z_vec[iz_low], digits=3))", 
           linewidth=2, xlabel="a", ylabel="σ(a,z)", title="Policy Function - Low Aversion (Different z)")
plot!(p2b, model_low_aversion.a_vec, σ_low[:, iz_mid], label="z=$(round(model_low_aversion.z_vec[iz_mid], digits=3))", 
      linewidth=2, linestyle=:dash)
plot!(p2b, model_low_aversion.a_vec, σ_low[:, iz_high], label="z=$(round(model_low_aversion.z_vec[iz_high], digits=3))", 
      linewidth=2, linestyle=:dot)
plot(p2b)

p2c = plot(model_high_aversion.a_vec, σ_high[:, iz_low], label="z=$(round(model_high_aversion.z_vec[iz_low], digits=3))", 
           linewidth=2, xlabel="a", ylabel="σ(a,z)", title="Policy Function - High Aversion (Different z)")
plot!(p2c, model_high_aversion.a_vec, σ_high[:, iz_mid], label="z=$(round(model_high_aversion.z_vec[iz_mid], digits=3))", 
      linewidth=2, linestyle=:dash)
plot!(p2c, model_high_aversion.a_vec, σ_high[:, iz_high], label="z=$(round(model_high_aversion.z_vec[iz_high], digits=3))", 
      linewidth=2, linestyle=:dot)
plot(p2c)

# Policy Functions - Zoomed in near a_min
p2b_zoom = plot(model_low_aversion.a_vec, σ_low[:, iz_low], label="z=$(round(model_low_aversion.z_vec[iz_low], digits=3))", 
                 linewidth=2, xlabel="a", ylabel="σ(a,z)", title="Policy Function near a_min (Low Aversion)")
plot!(p2b_zoom, model_low_aversion.a_vec, σ_low[:, iz_mid], label="z=$(round(model_low_aversion.z_vec[iz_mid], digits=3))", 
      linewidth=2, linestyle=:dash)
plot!(p2b_zoom, model_low_aversion.a_vec, σ_low[:, iz_high], label="z=$(round(model_low_aversion.z_vec[iz_high], digits=3))", 
      linewidth=2, linestyle=:dot)
xlims!(p2b_zoom, model_low_aversion.a_min, zoom_range_low)

p2c_zoom = plot(model_high_aversion.a_vec, σ_high[:, iz_low], label="z=$(round(model_high_aversion.z_vec[iz_low], digits=3))", 
                 linewidth=2, xlabel="a", ylabel="σ(a,z)", title="Policy Function near a_min (High Aversion)")
plot!(p2c_zoom, model_high_aversion.a_vec, σ_high[:, iz_mid], label="z=$(round(model_high_aversion.z_vec[iz_mid], digits=3))", 
      linewidth=2, linestyle=:dash)
plot!(p2c_zoom, model_high_aversion.a_vec, σ_high[:, iz_high], label="z=$(round(model_high_aversion.z_vec[iz_high], digits=3))", 
      linewidth=2, linestyle=:dot)
xlims!(p2c_zoom, model_high_aversion.a_min, zoom_range_high)


consumption_low_aversion = get_consumption_matrix(model_low_aversion, σ_low)
consumption_high_aversion = get_consumption_matrix(model_high_aversion, σ_high)

# Consumption function plots - Low aversion (full range)
p3b = plot(model_low_aversion.a_vec, consumption_low_aversion[:, iz_low], label="z=$(round(model_low_aversion.z_vec[iz_low], digits=3))", 
           linewidth=2, xlabel="a", ylabel="c(a,z)", title="Consumption Function - Low Aversion (Different z)")
plot!(p3b, model_low_aversion.a_vec, consumption_low_aversion[:, iz_mid], label="z=$(round(model_low_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p3b, model_low_aversion.a_vec, consumption_low_aversion[:, iz_high], label="z=$(round(model_low_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
display(p3b)

# Consumption function - zoomed near a_min (low aversion)
p3b_zoom = plot(model_low_aversion.a_vec, consumption_low_aversion[:, iz_low], label="z=$(round(model_low_aversion.z_vec[iz_low], digits=3))", 
                 linewidth=2, xlabel="a", ylabel="c(a,z)", title="Consumption near a_min (Low Aversion)")
plot!(p3b_zoom, model_low_aversion.a_vec, consumption_low_aversion[:, iz_mid], label="z=$(round(model_low_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p3b_zoom, model_low_aversion.a_vec, consumption_low_aversion[:, iz_high], label="z=$(round(model_low_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
xlims!(p3b_zoom, model_low_aversion.a_min, zoom_range_low)
display(p3b_zoom)

# Consumption function plots - High aversion (full range)
p3c = plot(model_high_aversion.a_vec, consumption_high_aversion[:, iz_low], label="z=$(round(model_high_aversion.z_vec[iz_low], digits=3))", 
           linewidth=2, xlabel="a", ylabel="c(a,z)", title="Consumption Function - High Aversion (Different z)")
plot!(p3c, model_high_aversion.a_vec, consumption_high_aversion[:, iz_mid], label="z=$(round(model_high_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p3c, model_high_aversion.a_vec, consumption_high_aversion[:, iz_high], label="z=$(round(model_high_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
display(p3c)

# Consumption function - zoomed near a_min (high aversion)
p3c_zoom = plot(model_high_aversion.a_vec, consumption_high_aversion[:, iz_low], label="z=$(round(model_high_aversion.z_vec[iz_low], digits=3))", 
                 linewidth=2, xlabel="a", ylabel="c(a,z)", title="Consumption near a_min (High Aversion)")
plot!(p3c_zoom, model_high_aversion.a_vec, consumption_high_aversion[:, iz_mid], label="z=$(round(model_high_aversion.z_vec[iz_mid], digits=3))", linewidth=2, linestyle=:dash)
plot!(p3c_zoom, model_high_aversion.a_vec, consumption_high_aversion[:, iz_high], label="z=$(round(model_high_aversion.z_vec[iz_high], digits=3))", linewidth=2, linestyle=:dot)
xlims!(p3c_zoom, model_high_aversion.a_min, zoom_range_high)
display(p3c_zoom)

#calculating euler residuals on a test grid Euler 
euler_residuals_low_aversion_grid, euler_residuals_low_aversion = euler_residuals(model_low_aversion, σ_low)
euler_residuals_high_aversion_grid, euler_residuals_high_aversion = euler_residuals(model_high_aversion, σ_high)

#plotting Euler residuals 
plot(euler_residuals_low_aversion, xlabel="a", ylabel="Euler Residuals", title="Euler Equation Residuals (Low Aversion)")

plot(euler_residuals_high_aversion, xlabel="a", ylabel="Euler Residuals", title="Euler Equation Residuals (High Aversion)")
