using Random, Distributions, StatsPlots, Statistics

Random.seed!(1234)

λ = 1.0
n = 1000

pois = Poisson(λ)

#draw_pois = rand(pois, n)

ranges = [5, 25, 100, 1000]
# needs corrections
for i in 1:4
    draw_pois = zeros(1000)
    for j in 1:ranges[i] #summing
        draw_pois += rand(pois, n)
    end
    pois_standarized = (draw_pois .- n .* λ)./(λ * n) #i dont remember if this standarization is ok?
end

#histogram looks okayish
histogram(pois_standarized;
	    bins=30,
	    normalize = :pdf,
	    label = "data",
	    xlabel = "pois_standarized",
	    ylabel = "Density",
	    title = "Histogram of pois_standarized")
    x = range(minimum(pois_standarized), maximum(pois_standarized), length=200)
    μ = mean(pois_standarized)
    σ = std(pois_standarized)
    
# kind of brute force below 
draw_pois = zeros(1000)
for j in 1:5 #summing
    draw_pois += rand(pois, n)
end
pois_standarized = (draw_pois .- n .* λ)./(λ * n) #i dont remember if this standarization is ok?
histogram(pois_standarized;
	    bins=30,
	    normalize = :pdf,
	    label = "data",
	    xlabel = "pois_standarized",
	    ylabel = "Density",
	    title = "Histogram of pois_standarized")
x = range(minimum(pois_standarized), maximum(pois_standarized), length=200)
μ = mean(pois_standarized)
σ = std(pois_standarized)

for j in 1:25 #summing
    draw_pois += rand(pois, n)
end
pois_standarized = (draw_pois .- n .* λ)./(λ * n) #i dont remember if this standarization is ok?
histogram(pois_standarized;
	    bins=30,
	    normalize = :pdf,
	    label = "data",
	    xlabel = "pois_standarized",
	    ylabel = "Density",
	    title = "Histogram of pois_standarized")
x = range(minimum(pois_standarized), maximum(pois_standarized), length=200)
μ = mean(pois_standarized)
σ = std(pois_standarized)

for j in 1:100 #summing
    draw_pois += rand(pois, n)
end
pois_standarized = (draw_pois .- n .* λ)./(λ * n) #i dont remember if this standarization is ok?
histogram(pois_standarized;
	    bins=30,
	    normalize = :pdf,
	    label = "data",
	    xlabel = "pois_standarized",
	    ylabel = "Density",
	    title = "Histogram of pois_standarized")
x = range(minimum(pois_standarized), maximum(pois_standarized), length=200)
μ = mean(pois_standarized)
σ = std(pois_standarized)

for j in 1:1000 #summing
    draw_pois += rand(pois, n)
end
pois_standarized = (draw_pois .- n .* λ)./(λ * n) #i dont remember if this standarization is ok?
histogram(pois_standarized;
	    bins=30,
	    normalize = :pdf,
	    label = "data",
	    xlabel = "pois_standarized",
	    ylabel = "Density",
	    title = "Histogram of pois_standarized")
x = range(minimum(pois_standarized), maximum(pois_standarized), length=200)
μ = mean(pois_standarized)
σ = std(pois_standarized)