using LaTeXStrings
using PrettyTables, Plots, LinearAlgebra, NLsolve, Roots

σ = 0.2
ω_1 = [1.0, 1.0]
ω_2 = [0.5, 1.5]

p_1
p_2 = 1.0

#consumption_people 1 - please check for mistakes
consumption_people1_good1(p_1, p_2, x, σ, ω_1, ω_2) = (( x^σ * p_1 ^ (-σ) ) * (p_1 * ω_1[1] + p_2 * ω_1[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))
consumption_people1_good2(p_1, p_2, x, σ, ω_1, ω_2) = (( (1 - x)^σ * p_2 ^ (-σ) ) * (p_1 * ω_1[1] + p_2 * ω_1[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))

#consumption_people 2
consumption_people2_good1(p_1, p_2, x, σ, ω_1, ω_2) = (( x^σ * p_1 ^ (-σ) ) * (p_1 * ω_2[1] + p_2 * ω_2[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))
consumption_people2_good2(p_1, p_2, x, σ, ω_1, ω_2) = (( (1 - x)^σ * p_2 ^ (-σ) ) * (p_1 * ω_2[1] + p_2 * ω_2[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))

