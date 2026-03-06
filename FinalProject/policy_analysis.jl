using Revise, QuantEcon, Statistics, Parameters, Interpolations, Optim, NumericalIntegration, LinearAlgebra, Plots

# Ładowanie modułu
include(joinpath(@__DIR__, "module.jl"))
using .Engine: InitializeModel, Solve_Model, GENERATE_K_GRID, OPTIMAL_H, OUTPUT, stationary_distribution, moments

println("Aktywne wątki: ", Threads.nthreads())

# KALIBRACJA
# Wstawić tutaj parametry (gamma, F, ps), które wyszły z estymacji SMM
calibrated_gamma = 0.05
calibrated_F = 0.01
calibrated_ps = 0.80

println("\n=== 1. ROZWIĄZYWANIE BAZOWEGO MODELU (tau = 0.0) ===")
model_base = InitializeModel(gamma=calibrated_gamma, F=calibrated_F, ps=calibrated_ps, tau=0.0, N_A=250, N_z=5)
V_init_base = zeros(model_base.N_A, model_base.N_z)
result_base = Solve_Model(model_base, V_init_base, 1e-6, 3000)
μ_base = stationary_distribution(model_base, result_base[3])
K_GRID = GENERATE_K_GRID(model_base, :polynomial, 5.0)

println("\n=== 2. ROZWIĄZYWANIE MODELU Z SUBSYDIUM (tau = 0.10) ===")
model_sub = InitializeModel(gamma=calibrated_gamma, F=calibrated_F, ps=calibrated_ps, tau=0.10, N_A=250, N_z=5)
V_init_sub = zeros(model_sub.N_A, model_sub.N_z)
result_sub = Solve_Model(model_sub, V_init_sub, 1e-6, 3000)
μ_sub = stationary_distribution(model_sub, result_sub[3])

# OBLICZANIE AGREGATÓW
function calc_aggregates(model, μ, result, K_GRID)
    K_agg = 0.0
    H_agg = 0.0
    Y_agg = 0.0
    E_K = 0.0; E_Z = 0.0; E_KZ = 0.0; Var_K = 0.0; Var_Z = 0.0
    
    for iz in 1:model.N_z
        z = model.z_vec[iz]
        for ik in 1:model.N_A
            k = K_GRID[ik]
            prob = μ[ik, iz]
            
            K_agg += prob * k
            H_agg += prob * OPTIMAL_H(model, k, z)
            Y_agg += prob * OUTPUT(model, k, z)
            
            E_K += prob * k
            E_Z += prob * z
            E_KZ += prob * k * z
        end
    end
    
    for iz in 1:model.N_z
        z = model.z_vec[iz]
        for ik in 1:model.N_A
            prob = μ[ik, iz]
            Var_K += prob * (K_GRID[ik] - E_K)^2
            Var_Z += prob * (z - E_Z)^2
        end
    end
    
    Cov_KZ = E_KZ - E_K * E_Z
    Corr_KZ = Cov_KZ / sqrt(Var_K * Var_Z)
    
    return K_agg, H_agg, Y_agg, Corr_KZ
end

K_b, H_b, Y_b, Corr_b = calc_aggregates(model_base, μ_base, result_base, K_GRID)
K_s, H_s, Y_s, Corr_s = calc_aggregates(model_sub, μ_sub, result_sub, K_GRID)

# KOSZT SUBSYDIUM
subsidy_cost = 0.0
i_sub = result_sub[2]
for iz in 1:model_sub.N_z
    for ik in 1:model_sub.N_A
        inv = i_sub[ik, iz]
        if inv > 0
            subsidy_cost += μ_sub[ik, iz] * model_sub.tau * inv
        end
    end
end
cost_to_Y = subsidy_cost / Y_s

# WYPISYWANIE WYNIKÓW
println("\n=== PORÓWNANIE WYNIKÓW (BAZOWY vs SUBSYDIUM) ===")
mom_b = moments(μ_base, result_base, K_GRID)
mom_s = moments(μ_sub, result_sub, K_GRID)

println("1. Momenty Inwestycyjne:")
println("   Średnia inwestycja (i/k): $(round(mom_b[1]*100, digits=2))% -> $(round(mom_s[1]*100, digits=2))%")
println("   Inaction Rate: $(round(mom_b[2]*100, digits=2))% -> $(round(mom_s[2]*100, digits=2))%")
println("   Positive Spikes: $(round(mom_b[4]*100, digits=2))% -> $(round(mom_s[4]*100, digits=2))%")

println("\n2. Agregaty Makroekonomiczne:")
println("   Kapitał (K): $(round(K_b, digits=3)) -> $(round(K_s, digits=3)) ($(round((K_s/K_b - 1)*100, digits=2))%)")
println("   Praca (H): $(round(H_b, digits=3)) -> $(round(H_s, digits=3)) ($(round((H_s/H_b - 1)*100, digits=2))%)")
println("   Produkcja (Y): $(round(Y_b, digits=3)) -> $(round(Y_s, digits=3)) ($(round((Y_s/Y_b - 1)*100, digits=2))%)")

println("\n3. Przekrój poprzeczny (Cross-sectional):")
println("   Korelacja Corr(k, z): $(round(Corr_b, digits=4)) -> $(round(Corr_s, digits=4))")

println("\n4. Koszt Subsydium:")
println("   Całkowity koszt: $(round(subsidy_cost, digits=4))")
println("   Procent produkcji (Cost/Y): $(round(cost_to_Y*100, digits=2))%")

# WYKRESY
z_mid = Int(floor(model_base.N_z / 2)) + 1

# 1. Policy functions: i*(k, z)
p1 = plot(K_GRID[1:100], result_base[2][1:100, z_mid], label="Bazowy (tau=0)", title="Inwestycje i(k, z_mid)", lw=2)
plot!(p1, K_GRID[1:100], result_sub[2][1:100, z_mid], label="Subsydium (tau=0.10)", lw=2, linestyle=:dash)
hline!(p1, [0], color=:black, alpha=0.3, label="")
xlabel!("Kapitał k")
ylabel!("Inwestycja i")

# 2. Rozkład marginalny kapitału
μ_K_b = sum(μ_base, dims=2)[:, 1]
μ_K_s = sum(μ_sub, dims=2)[:, 1]
p2 = plot(K_GRID[1:100], μ_K_b[1:100], label="Bazowy (tau=0)", title="Rozkład stacjonarny K", lw=2)
plot!(p2, K_GRID[1:100], μ_K_s[1:100], label="Subsydium (tau=0.10)", lw=2, linestyle=:dash)
xlabel!("Kapitał k")
ylabel!("Prawdopodobieństwo")

final_plot = plot(p1, p2, layout=(2,1), size=(800, 800), margin=5Plots.mm)
display(final_plot)
savefig(final_plot, "policy_analysis_results.png")
println("\nWykresy zapisane do 'policy_analysis_results.png'")