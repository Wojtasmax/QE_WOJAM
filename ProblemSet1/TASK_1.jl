using Random, Distributions, StatsPlots, Statistics

Random.seed!(1234)

λ = 1.0
n = 1000

pois = Poisson(λ)

#draw_pois = rand(pois, n)

ranges = [5, 25, 100, 1000]
# needs corrections - draws only one plot - why?
for N in ranges
    draw_pois = zeros(1000)
    for j in 1:N #summing
        draw_pois += rand(pois, n)
    end
        h = histogram(pois_standarized;
        bins=20,
        normalize=:pdf,
        label="data",
        xlabel="pois_standarized",
        ylabel="Density",
        title="Histogram for N=$(N)")

    x = range(minimum(pois_standarized), maximum(pois_standarized), length=200)
    μ = mean(pois_standarized)
    σ = std(pois_standarized)
    plot!(h, x, pdf.(Normal(μ, σ), x); lw=2, color=:red, label="Normal(μ,σ)")

    display(h)    
end

pois_standarized
histogram(pois_standarized)



#here i tried to store values of each iteration  in one matrix but concatenation crashes
mat = zeros(1000,1)
typeof(mat)
draw_pois = zeros(1000)
for i in 1:4
    for j in 1:ranges[i] #summing
        draw_pois += rand(pois, n)
    end
    pois_standarized = (draw_pois .- n .* λ)./(λ * n) #i dont remember if this standarization is ok?
    hcat(mat, pois_standarized)
end
print(mat)

    
# kind of brute force below for ranges - it works


draw_pois = zeros(n)
for j in 1:n #summing
    draw_pois += rand(pois, n)
end

pois_standarized = (draw_pois./5 .- λ)./sqrt(λ / 5) #i dont remember if this standarization is ok?
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
plot!(x, pdf.(Normal(μ, σ), x), label="Normal(μ, σ)")



draw_pois = zeros(n)
for j in 1:25 #summing
    draw_pois += rand(pois, n)
end

pois_standarized = (draw_pois./25 .- λ)./sqrt(λ / 25) #i dont remember if this standarization is ok?
histogram(pois_standarized;
	    bins=20,
	    normalize = :pdf,
	    label = "data",
	    xlabel = "pois_standarized",
	    ylabel = "Density",
	    title = "Histogram of pois_standarized")
x = range(minimum(pois_standarized), maximum(pois_standarized), length=200)
μ = mean(pois_standarized)
σ = std(pois_standarized)
plot!(x, pdf.(Normal(μ, σ), x), label="Normal(μ, σ)")

draw_pois = zeros(n)
for j in 1:100  #summing
    draw_pois += rand(pois, n)
end

pois_standarized = (draw_pois./100 .- λ)./sqrt(λ / 100) #i dont remember if this standarization is ok?
histogram(pois_standarized;
	    bins=20,
	    normalize = :pdf,
	    label = "data",
	    xlabel = "pois_standarized",
	    ylabel = "Density",
	    title = "Histogram of pois_standarized")
x = range(minimum(pois_standarized), maximum(pois_standarized), length=200)
μ = mean(pois_standarized)
σ = std(pois_standarized)
plot!(x, pdf.(Normal(μ, σ), x), label="Normal(μ, σ)")


draw_pois = zeros(n)
for j in 1:1000  #summing
    draw_pois += rand(pois, n)
end

pois_standarized = (draw_pois./1000 .- λ)./sqrt(λ / 1000) #i dont remember if this standarization is ok?
histogram(pois_standarized;
	    bins=20,
	    normalize = :pdf,
	    label = "data",
	    xlabel = "pois_standarized",
	    ylabel = "Density",
	    title = "Histogram of pois_standarized")
x = range(minimum(pois_standarized), maximum(pois_standarized), length=200)
μ = mean(pois_standarized)
σ = std(pois_standarized)
plot!(x, pdf.(Normal(μ, σ), x), label="Normal(μ, σ)")


