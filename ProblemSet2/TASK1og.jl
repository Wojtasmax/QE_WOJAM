using QuadGK, Distributions, Roots, Plots, NLsolve
function foc_integral(ω, W, Rf, γ, μ, σ)
try
    return quadgk(r -> (r-Rf)*(W*(ω*r+(1-ω)*Rf))^(-γ)*pdf(LogNormal(μ, σ), r), 0, Inf)[1]
catch e
    return NaN
end
end
foc_integral(10, 1, 1.02, 2.0, 0.05, 0.1)
# test it for γ=0
# first naive check:
function expected_result(μ, σ, Rf)
    return exp(μ+0.5*σ^2)-Rf
end
function are_they_equal()
    for i in 1:20
        ω, W, Rf, μ, σ=rand(), 200*rand(), 3*rand()+1,rand(), 2*rand()
        if !isapprox(foc_integral(ω, W, Rf, 0.0, μ, σ), expected_result(μ,σ,Rf); atol=1e-8, rtol=1e-6)
                return false
                break
        end
    end
    return true
end
print(are_they_equal() ? "The functions are approximately equal\n" : "The functions differ substantially\n")
# 1.3
function optimal_portfolio(W, Rf, γ, μ, σ)
    f(ω) = foc_integral(ω, W, Rf, γ, μ, σ)
    ω_star=find_zero(f, 0.5, Order1())
    return ω_star
end
function optimal_portfolio(W, Rf, γ, μ, σ)
    if abs(γ)<1
        return 1.
    end
    function f!(F, x)
        ω = x[1]
        F[1] = foc_integral(ω, W, Rf, γ, μ, σ)
    end
    initial_guess = [0.5]
    return nlsolve(f!, initial_guess).zero[1]
end
# test it for some parameters
optimal_portfolio(1, 1.02, 3.0, 0.05, 0.1)
# 1.5
const_optimal_portfolio(γ) = optimal_portfolio(1, 1.02, γ, 0.05, 0.1)
plot(const_optimal_portfolio, .1, 10.0, xlabel="Risk Aversion γ", ylabel="Optimal Share in Risky Asset ω*",
    title="Optimal Portfolio Share vs Risk Aversion", legend=false)
# the chart cuts at around 3. thats when the optimal share is above 1
const_optimal_portfolio(2)