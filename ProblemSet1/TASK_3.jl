using LaTeXStrings
using PrettyTables, Plots, LinearAlgebra, NLsolve, Roots

σ = 0.2
ω_1 = [1.0, 1.0]
ω_2 = [0.5, 1.5]

function demand_good1(p1, x, σ, ω)
    num = (x^σ) * p1^(-σ) * (p1*ω[1] + ω[2])
    den = x^σ * p1^(1 - σ) + (1 - x)^σ
    return num / den
end

function solve(share, elasticity, wealth1, wealth2) 

    equation(p1) =  demand_good1(p1, share, elasticity, wealth1) + demand_good1(p1, share, elasticity, wealth2) - (wealth1[1] + wealth2[1])
    #Solving for logaritm, because we deal with negative numbers - z can be from -Inf to +Inf, while p must be greater than 0.
    equation_z(z) = equation(exp(z))
    
    try
        zstar = find_zero(equation_z, (-100.0, 100.0), Bisection())
        return exp(zstar)
    catch e
        # if solution is extreme
        return NaN
    end

    return exp(zstar)
end 

#Choosing the range of x
shares = range(0.01, 0.99; length = 1000)

#Solving for solutions for x in range (0,1)
solutions_price = [
    solve(x, σ, ω_1, ω_2)
    for x in shares
]

#Ploting the prices of good 1 
h1 = plot(shares, solutions_price, label = "Price of good 1, σ = 0.2")

solutions_consumption = [
    demand_good1(solve(x, σ, ω_1, ω_2), x, σ, ω_1)
    for x in shares
]

#Ploting the consumption of good 1 by both people with old elasticity
plot(shares, solutions_consumption, label = "Consumption of good 1 by person 1, σ = 0.2")
h2 = plot!(shares, ω_1[1]+ω_2[1].-solutions_consumption, label = "Consumption of good 1 by person 2, σ = 0.2")

#Changing the elasticity to 5.0
σ_2 = 5.0

solutions_price_2 = [
    solve(x, σ_2, ω_1, ω_2)
    for x in shares
]

#comparing prices - needs domain improvement, and more estethic taste
plot(shares, solutions_price, maximum =1.1, label = "σ = 0.2", xlabel ="share", ylabel = "price")
h3 = plot!(shares, solutions_price_2, maximum = 1.1, label = "σ = 5.0 ", xlabel ="share", ylabel = "price")

solutions_consumption2 = [
    demand_good1(solve(x, σ_2, ω_1, ω_2), x, σ, ω_1)
    for x in shares
]
  
#Ploting the consumption of good 1 by both people with new elasticity
plot(shares, solutions_consumption2, label = "Consumption of good 1 by person 1 with σ = 5")
h4 = plot!(shares, ω_1[1]+ω_2[1].-solutions_consumption2, label = "Consumption of good 1 by person 2 with σ = 5")




#The higher the elasticity, the more substitute the goods become. Therefore, given the starting parameters, the growing elasticity would mean that they would both prefer to have just one good, and so price would not be very dependent on x. 
#With lower elasticity they will greatly benefit from having two goods and therefore the equilibrium price will be more dependent on x (because it will greatly affect what proportions they desire)
#Obviously for the situations of x=0 and x=1 all the functions yield the same result, but that is generally outside the CES normal aplications.