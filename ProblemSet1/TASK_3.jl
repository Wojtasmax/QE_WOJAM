using LaTeXStrings
using PrettyTables, Plots, LinearAlgebra, NLsolve, Roots

σ = 0.2
ω_1 = [1.0, 1.0]
ω_2 = [0.5, 1.5]
p_2 = 1.0 #For simplicity, we assume good 2 is numeraire

consumption_people1_good1(p_1; p_2 , x, σ, ω_1, ω_2) = (( x^σ * p_1 ^ (-σ) ) * (p_1 * ω_1[1] + p_2 * ω_1[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))
consumption_people1_good2(p_1; p_2, x, σ, ω_1, ω_2) = (( (1 - x)^σ * p_2 ^ (-σ) ) * (p_1 * ω_1[1] + p_2 * ω_1[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))

#Do we have to define it twice, and later on 
consumption_people2_good1(p_1; p_2, x, σ, ω_1, ω_2) = (( x^σ * p_1 ^ (-σ) ) * (p_1 * ω_2[1] + p_2 * ω_2[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))
consumption_people2_good2(p_1; p_2, x, σ, ω_1, ω_2) = (( (1 - x)^σ * p_2 ^ (-σ) ) * (p_1 * ω_2[1] + p_2 * ω_2[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))

function solve(price, share, elasticity, wealth1, wealth2) 
    
    consumption_people1_good1(p_1; p_2 , x, σ, ω_1, ω_2) = (( x^σ * p_1 ^ (-σ) ) * (p_1 * ω_1[1] + p_2 * ω_1[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))
    consumption_people1_good2(p_1; p_2, x, σ, ω_1, ω_2) = (( (1 - x)^σ * p_2 ^ (-σ) ) * (p_1 * ω_1[1] + p_2 * ω_1[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))

    #consumption_people 2
    consumption_people2_good1(p_1; p_2, x, σ, ω_1, ω_2) = (( x^σ * p_1 ^ (-σ) ) * (p_1 * ω_2[1] + p_2 * ω_2[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))
    consumption_people2_good2(p_1; p_2, x, σ, ω_1, ω_2) = (( (1 - x)^σ * p_2 ^ (-σ) ) * (p_1 * ω_2[1] + p_2 * ω_2[2]))/(x^σ * p_1^(1 - σ) + (1 - x)^σ * p_2 ^(1 - σ))

    equation(p_1) =  consumption_people1_good1(p_1; p_2 = price, x = share,σ = elasticity, ω_1 = wealth1, ω_2 = wealth2) + consumption_people2_good1(p_1; p_2 = price, x = share,σ = elasticity, ω_1 = wealth1, ω_2 = wealth2 ) - (wealth1[1] + wealth2[1])
    equation_z(z) = equation(exp(z))
    zstar = find_zero(equation_z, 0.0; autodiff = :forward)
    return exp(zstar)
end 

#Choosing the range of x
shares = range(0.01, 0.99; length = 1000)

#Solving for solutions for x in range (0,1)
solutions_price = [
    solve(p_2, x, σ, ω_1, ω_2)
    for x in shares
]

#Ploting the solutions
plot(shares, solutions_price)

solutions_consumption = [
    consumption_people1_good1(solve(p_2, x, σ, ω_1, ω_2); p_2, x, σ, ω_1, ω_2)
    for x in shares
]

plot(shares, solutions_consumption, label = "Consumption of good 1 by person 1")
plot!(shares, ω_1[1]+ω_2[1].-solutions_consumption, label = "Consumption of good 1 by person 2")

#Changing the elasticity
σ_2 = 1.52

solutions_price_2 = [
    solve(p_2, x, σ_2, ω_1, ω_2)
    for x in shares
]

#Ploting the solutions
plot(shares, solutions_price)
plot!(shares, solutions_price_2)
