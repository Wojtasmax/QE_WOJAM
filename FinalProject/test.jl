

###TO JEST TYLKO DO SPRAWDZENIA CZY TO SIE WOGOLE LICZY!!!!!!!!!!!!!!!!!!!!#############











using QuantEcon, Statistics, Parameters, Interpolations, Optim
# Załaduj moduł Engine z tego samego katalogu
include(joinpath(@__DIR__, "module.jl"))
using .Engine

# Inicjalizacja modelu
println("Inicjalizacja modelu...")
model_base = ProjectParams()
model = ProjectParams(model_base)

println("Parametry modelu:")
println("α = ", model.α)
println("β = ", model.β)
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

# Rozwiązanie modelu
result = Solve_Model(model, V_init, 1e-4, 10000)

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

# 1. Przygotowanie danych (ręczne wyliczenie gridu dla pewności)
K_GRID = model.k_min .+ (model.k_max - model.k_min) .* (model.ω .^ 5.0)

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