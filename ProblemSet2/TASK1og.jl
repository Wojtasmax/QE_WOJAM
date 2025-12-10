using QuadGK, Distributions, Roots, Plots
function foc_integral(ω, W, Rf, γ, μ, σ)
    return quadgk(r -> (r-Rf)*(W*(ω*r+(1-ω)*Rf))^(-γ)*pdf(LogNormal(μ, σ), r), 0, Inf)[1]
end
# test it for γ=0
# first naive check:
for i in 1:5
ω, W, Rf, μ, σ=rand(), 200*rand(), 3*rand()+1,rand(), 2*rand()
print("Approximate error=",foc_integral(ω, W, Rf, 0.0, μ, σ)-(exp(μ+0.5*σ^2)-Rf), "\n")
end
# now build in julia function
println(isapprox(foc_integral(ω, W, Rf, 0.0, μ, σ), exp(μ+0.5*σ^2)-Rf, atol=1e-5) ? "The functions are almost equal" : "The functions differ substantially")
# 1.3
function optimal_portfolio(W, Rf, γ, μ, σ)
    f(ω) = foc_integral(ω, W, Rf, γ, μ, σ)
    ω_star=find_zero(f, 0.5, Order1())
    return ω_star
end
# test it for some parameters
optimal_portfolio(1, 1.02, 3.0, 0.05, 0.1)
# 1.5
const_optimal_portfolio(γ) = optimal_portfolio(1, 1.02, γ, 0.05, 0.1)
plot(const_optimal_portfolio, .1, 10.0, xlabel="Risk Aversion γ", ylabel="Optimal Share in Risky Asset ω*",
    title="Optimal Portfolio Share vs Risk Aversion", legend=false)
# the chart cuts at around 3. thats when the optimal share is above 1