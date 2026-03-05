### TO JEST TYLKO DO SPRAWDZENIA CZY TO SIE WOGOLE LICZY!!!!!!!!!!!!!!!!!!!!#############

using Revise, QuantEcon, Statistics, Parameters, Interpolations, Optim, NumericalIntegration, LinearAlgebra

# Załaduj moduł Engine z tego samego katalogu
include(joinpath(@__DIR__, "module.jl"))
using .Engine: InitializeModel, Solve_Model, VFI, GENERATE_K_GRID, ADJ_COST, IRR, OPTIMAL_H, OPERATING_PROFIT, get_transition_matrix_young, stationary_distribution, moments, objective_function_SMM

# Sprawdzenie wątków (Upewnij się, że odpaliłeś Julię przez julia -t auto)
println("Aktywne wątki: ", Threads.nthreads())

# Inicjalizacja modelu
println("Inicjalizacja modelu...")
model = InitializeModel()  # ZMIANA: używamy InitializeModel zamiast InitilizeModel

println("Parametry modelu:")
println("α = ", model.alpha)
println("r = ", model.r)
println("β = ", 1.0 / (1.0 + model.r))  # ZMIANA: obliczamy beta na miejscu
println("N_z = ", model.N_z)
println("N_A = ", model.N_A)
println()

# Inicjalizacja value function
N_A = model.N_A
N_z = model.N_z
V_init = zeros(N_A, N_z)

println("Rozpoczynam rozwiązywanie modelu...")
println("Wymiary V_init: ", size(V_init))
println()

# Rozwiązanie modelu - TU ODPALA SIĘ MULTITHREADING Z TWOJEGO MODULE.JL
result = Solve_Model(model, V_init, 1e-12, 10000)

println()
println("Model rozwiązany!")
println("Maksymalna wartość funkcji wartości: ", maximum(result[1]))
println("Minimalna wartość funkcji wartości: ", minimum(result[1]))
println()

println("Statystyki polityki inwestycyjnej:")
println("Średnia inwestycja: ", mean(result[2]))
println("Max inwestycja: ", maximum(result[2]))
println("Min inwestycja: ", minimum(result[2]))
println()

println("Statystyki operating profit:")
println("Średni zysk operacyjny: ", mean(result[6]))
println("Max zysk operacyjny: ", maximum(result[6]))
println("Min zysk operacyjny: ", minimum(result[6]))
println()

println("GOTOWE! Wyniki zapisane w zmiennej 'result'")

using Plots

# 1. Przygotowanie danych - ZMIANA: generujemy K_GRID przy użyciu funkcji
K_GRID = GENERATE_K_GRID(model, :polynomial, 5.0)

z_low  = 1
z_mid  = Int(floor(model.N_z / 2)) + 1
z_high = model.N_z

# --- GENEROWANIE WYKRESÓW ---

# Wykres 1: Funkcja Wartości
p1 = plot(K_GRID, result[1][:, z_low], label="Niska produktywność", title="Funkcja Wartości V(k, z)", lw=2)
plot!(p1, K_GRID, result[1][:, z_mid], label="Średnia produktywność", lw=2)
plot!(p1, K_GRID, result[1][:, z_high], label="Wysoka produktywność", lw=2)
ylabel!("Wartość firmy V")

# Wykres 2: Decyzja o kapitale na jutro
p2 = plot(K_GRID, result[3][:, z_low], label="z_low", title="Polityka Kapitałowa k'(k, z)", lw=2)
plot!(p2, K_GRID, result[3][:, z_mid], label="z_mid", lw=2)
plot!(p2, K_GRID, result[3][:, z_high], label="z_high", lw=2)
plot!(p2, K_GRID, K_GRID, label="k' = k (Steady State)", linestyle=:dash, color=:black, alpha=0.5)
ylabel!("Kapitał jutro k'")

# Wykres 3: Inwestycje
p3 = plot(K_GRID, result[2][:, z_low], label="z_low", title="Inwestycje i(k, z)", lw=2)
plot!(p3, K_GRID, result[2][:, z_mid], label="z_mid", lw=2)
plot!(p3, K_GRID, result[2][:, z_high], label="z_high", lw=2)
hline!(p3, [0], label="", color=:black, alpha=0.3)
xlabel!("Kapitał dzisiaj k")
ylabel!("Inwestycja i")

# Złożenie i zapis
final_plot = plot(p1, p2, p3, layout=(3,1), size=(800, 1200), margin=5Plots.mm)
display(final_plot)
savefig(final_plot, "wyniki_modelu_firmy.png")

println("\n--- INTERPRETACJA WYNIKÓW ---")
println("1. Funkcja Wartości (V): Pokazuje całkowitą zdyskontowaną wartość firmy. Jest rosnąca i wklęsła względem kapitału, co oznacza malejące krańcowe korzyści z kapitału.")
println("2. Polityka Kapitałowa (k'): Pokazuje, ile kapitału firma chce mieć w następnym okresie. Punkt przecięcia z linią 45 stopni (przerywana) to stan stacjonarny - tam kapitał przestaje rosnąć/maleć.")
println("3. Inwestycje (i): Pokazuje bieżące zakupy/sprzedaż kapitału. Zwróć uwagę na płaskie odcinki (i=0) - to regiony, w których koszty dostosowania (F, gamma) sprawiają, że firmie nie opłaca się zmieniać kapitału.")

#### Stationary distribution

# Stationary distribution calculation
μ = stationary_distribution(model, result[3])

println("\nRozkład stacjonarny μ obliczony.")
println("Suma μ: ", sum(μ), " (powinno być ≈ 1.0)")

# Marginal distributions
μ_K = sum(μ, dims=2)[:, 1]  # Marginal distribution over capital
μ_z = sum(μ, dims=1)[1, :]  # Marginal distribution over productivity

# Investment rate
investment_rate = result[2] ./ K_GRID

# 5. Stationary distribution: Distribution over (k, z) and histogram of investment rates i/k. 
# This is before calculating method of moments, so we can see the raw distribution of investment rates without summarizing it into a single number.

p_dist1 = histogram(K_GRID, weights=μ_K, title="Marginal Distribution over Capital", 
                    xlabel="Capital k", ylabel="Probability", bins=30)

p_dist2 = bar(model.z_vec, μ_z, title="Marginal Distribution over Productivity", 
              xlabel="Productivity z", ylabel="Probability", legend=false)

p_dist3 = histogram(investment_rate[:], weights=μ[:], 
                    title="Distribution of Investment Rates", 
                    xlabel="Investment Rate i/k", ylabel="Frequency", 
                    bins=50, alpha=0.7, legend=false)

dist_plot = plot(p_dist1, p_dist2, p_dist3, layout=(3,1), size=(800, 1000))
display(dist_plot)
savefig(dist_plot, "stationary_distributions.png")

# Momenty
model_moments = moments(μ, result, K_GRID)

theta_init = [0.05, 0.01, 0.80]
opt_result = optimize(objective_function_SMM, theta_init, NelderMead(), Optim.Options(iterations=30, show_trace=true))

println("Skonfigurowane paramtery: ", opt_result)

println("\n=== ANALIZA ZAKOŃCZONA ===")
println("Wyniki zapisane w zmiennych:")
println("  - result: pełne wyniki VFI")
println("  - μ: rozkład stacjonarny")
println("  - model_moments: momenty")