using Random, Distributions, StatsPlots, Statistics

#setting the seed for reproducibility
Random.seed!(1234)

# setting the parameters
λ = 1.0
N = 1000

# defining Poisson distribution
pois = Poisson(λ)

# preparing array to store plots
plots = Vector{Plots.Plot}(undef, 4)

# simulating and plotting for different n
ranges = [5, 25, 100, 1000]

#loop over different n
for n in 1:4
	
	#preparing vector to store Poisson draws
	draw_pois = zeros(N)

	#summing n Poisson draws
    for j in 1:ranges[n] 
        draw_pois += rand(pois, N)
    end
    
	#standardizing mean with central limit theorem
	pois_standardized = (draw_pois./ranges[n] .-  λ)./(λ / sqrt(ranges[n])) 

	#plotting histogram
	h = histogram(pois_standardized;
        bins=25,
        normalize=:pdf,
        label="data",
        xlabel="Poisson standardized",
        ylabel="Density",
        title="Histogram for n=$(ranges[n])")

	#comparing with normal distribution
    x = range(minimum(pois_standardized), maximum(pois_standardized), length=200)
    μ = mean(pois_standardized)
    σ = std(pois_standardized)
    plot!(h, x, pdf.(Normal(μ, σ), x); lw=2, color=:red, label="Normal(μ,σ)")
    plots[n] = h   
end

#displaying all plots together
plot(plots..., layout=(2,2), size=(800,600)) 

#saving the figure
savefig("poisson_central_limit_theorem.png")