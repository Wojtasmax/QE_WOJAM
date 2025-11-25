using Random, Distributions, StatsPlots, Statistics

Random.seed!(1234)

λ = 1.0
N = 1000

pois = Poisson(λ)

plots = Vector{Plots.Plot}(undef, 4)

ranges = [5, 25, 100, 1000]

for n in 1:4
    draw_pois = zeros(N)
    for j in 1:ranges[n] #summing
        draw_pois += rand(pois, N)
    end
    pois_standardized = (draw_pois./ranges[n] .-  λ)./(λ / sqrt(ranges[n])) #needs to be verified
    h = histogram(pois_standardized;
        bins=25,
        normalize=:pdf,
        label="data",
        xlabel="Poisson standardized",
        ylabel="Density",
        title="Histogram for n=$(ranges[n])")


    x = range(minimum(pois_standardized), maximum(pois_standardized), length=200)
    μ = mean(pois_standardized)
    σ = std(pois_standardized)
    plot!(h, x, pdf.(Normal(μ, σ), x); lw=2, color=:red, label="Normal(μ,σ)")
    plots[n] = h   
end

plot(plots..., layout=(2,2), size=(800,600))

savefig("poisson_central_limit_theorem.png")