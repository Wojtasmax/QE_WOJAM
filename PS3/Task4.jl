using Parameters, Plots, LinearAlgebra, QuantEcon, Random, Distributions
Random.seed!(1234)

# =====================================================
# Setup
# =====================================================

Setup = @with_kw (
    # Education grid
    e_grid_min = 0.0,           # smallest grid point
    e_grid_max = 1.0,           # highest grid point
    e_grid_size = 100,          # grid size
    
    # Wage process
    w_grid_size = 7,            # grid size to discretize the AR(1) process 
    σ_ϵ = 0.15,                 # stdev of the AR(1) process for wages
    ϱ = 0.98,                   # persistence of the AR(1) process 
    w = 1.0,                    # avarage wage level in the economy
    
    # Human capital grid
    δ = 0.1,                    # depreciation rate
    h_grid_min = 0.0,           # smallest grid point
    h_grid_max = 2.5,           # we set the h_max to 10
    h_grid_size = 1000,         # grid size
)

Agent = @with_kw (
    β = 0.95,                   
    σ = 1.5,
    Ψ = 0.5,
    α = 0.1,
    f̂ = 5.0,
    u = σ == 1 ? log : c -> (c^(1-σ)-1)/(1-σ)
)

setup = Setup()
ag = Agent()

# =====================================================
# Grids
# =====================================================

e_grid = collect(range(setup.e_grid_min, setup.e_grid_max, length=setup.e_grid_size))
h_grid = collect(range(setup.h_grid_min, setup.h_grid_max, length=setup.h_grid_size))

# Human capital production
f = h -> min(h^ag.α + 0.1, ag.f̂)

# =====================================================
# Wage process (Rouwenhorst)
# =====================================================

mc = rouwenhorst(setup.w_grid_size, setup.ϱ, setup.σ_ϵ)
w_grid = setup.w .* exp.(mc.state_values)
P = mc.p

# =====================================================
# Allocate arrays
# =====================================================

V  = zeros(setup.h_grid_size, setup.w_grid_size)
Vn = similar(V)

e_policy = zeros(setup.h_grid_size, setup.w_grid_size)

# =====================================================
# Value Function Iteration
# =====================================================

tol = 1e-6
maxiter = 1000

for it in 1:maxiter
    for ih in 1:setup.h_grid_size
        h = h_grid[ih]

        for iw in 1:setup.w_grid_size
            w = w_grid[iw]

            v_best = -Inf
            e_best = 0.0

            for e in e_grid

                # Income
                y = w * f(h) * (1 - e)
                if y <= 0
                    continue
                end

                # Utility
                u_now = ag.u(y) - ag.Ψ * e

                # Human capital transition
                h_next = clamp(h + e - setup.δ,
                               setup.h_grid_min,
                               setup.h_grid_max)

                ih_next = searchsortedfirst(h_grid, h_next)
                ih_next = clamp(ih_next, 1, setup.h_grid_size)

                # Expected continuation value
                EV = dot(P[iw, :], V[ih_next, :])

                v = u_now + ag.β * EV

                if v > v_best
                    v_best = v
                    e_best = e
                end
            end

            Vn[ih, iw] = v_best
            e_policy[ih, iw] = e_best
        end
    end

    if maximum(abs.(Vn .- V)) < tol
        println("Converged in iteration $it")
        break
    end

    V .= Vn
end

# =====================================================
# Plots
# =====================================================

plot1 = plot(h_grid, e_policy[:,1], label="Low wage", lw=2)
plot!(h_grid, e_policy[:,4], label="Mid wage", lw=2)
plot!(h_grid, e_policy[:,7], label="High wage", lw=2)
xlabel!("Human capital h")
ylabel!("Education effort e")
title!("Education policy")


plot2 = plot(h_grid, V[:,1], label="Low wage", lw=2)
plot!(h_grid, V[:,4], label="Mid wage", lw=2)
plot!(h_grid, V[:,7], label="High wage", lw=2)
xlabel!("Human capital h")
ylabel!("Value V(h,w)")
title!("Value function for different wage states")


y_low  = zeros(setup.h_grid_size)
y_mid  = zeros(setup.h_grid_size)
y_high = zeros(setup.h_grid_size)

for ih in 1:setup.h_grid_size
    h = h_grid[ih]

    y_low[ih]  = w_grid[1] * f(h) * (1 - e_policy[ih, 1])
    y_mid[ih]  = w_grid[4] * f(h) * (1 - e_policy[ih, 4])
    y_high[ih] = w_grid[7] * f(h) * (1 - e_policy[ih, 7])
end

plot3 = plot(h_grid, y_low, label="Low wage", lw=2)
plot!(h_grid, y_mid, label="Mid wage", lw=2)
plot!(h_grid, y_high, label="High wage", lw=2)
xlabel!("Human capital h")
ylabel!("Income y")
title!("Income as a function of human capital")

# =====================================================
# Simulation parameters
# =====================================================

function simulate_path(T, h0, w0_index, h_grid, w_grid, e_policy, P, setup, f)
    h = zeros(T)
    w = zeros(T)
    e = zeros(T)
    y = zeros(T)
    w_index = zeros(Int, T)

    h[1] = h0
    w_index[1] = w0_index
    w[1] = w_grid[w0_index]

    for t in 1:T-1
        ih = searchsortedfirst(h_grid, h[t])
        ih = clamp(ih, 1, setup.h_grid_size)

        iw = w_index[t]

        e[t] = e_policy[ih, iw]
        y[t] = w[t] * f(h[t]) * (1 - e[t])

        h[t+1] = clamp(h[t] + e[t] - setup.δ,
                        setup.h_grid_min,
                        setup.h_grid_max)

        w_index[t+1] = rand(Categorical(P[iw, :]))
        w[t+1] = w_grid[w_index[t+1]]
    end

    # last period
    ih = searchsortedfirst(h_grid, h[T])
    ih = clamp(ih, 1, setup.h_grid_size)
    iw = w_index[T]
    e[T] = e_policy[ih, iw]
    y[T] = w[T] * f(h[T]) * (1 - e[T])

    return h, w, e, y
end


T = 1100;
burnin = 100;
obs = 101:200
keep = 201:1100

Nsim = 5

h_sims = Vector{Vector{Float64}}(undef, Nsim)
w_sims = Vector{Vector{Float64}}(undef, Nsim)
e_sims = Vector{Vector{Float64}}(undef, Nsim)
y_sims = Vector{Vector{Float64}}(undef, Nsim)

w_keep = Vector{Vector{Float64}}(undef, Nsim)
e_keep = Vector{Vector{Float64}}(undef, Nsim)
y_keep = Vector{Vector{Float64}}(undef, Nsim)

y_all = Float64[]
e_all = Float64[]
w_all = Float64[]

for n in 1:Nsim
    h, w, e, y = simulate_path(
        T,
        1.0,              # h0
        4,                # w4
        h_grid,
        w_grid,
        e_policy,
        P,
        setup,
        f
    )

    h_sims[n] = h[obs]
    w_sims[n] = w[obs]
    e_sims[n] = e[obs]
    y_sims[n] = y[obs]

    w_keep[n] = w[keep]
    e_keep[n] = e[keep]
    y_keep[n] = y[keep]

end

for n in 1:Nsim
    append!(y_all, y_keep[n][1:900])
    append!(e_all, e_keep[n][1:900])
    append!(w_all, w_keep[n][1:900])
end

plot4 = plot()
for n in 1:Nsim
    plot!(h_sims[n], label="Sim $n", lw=2)
end
xlabel!("Time (101–200)")
ylabel!("Human capital h")
title!("Human capital paths")

plot5 = plot()
for n in 1:Nsim
    plot!(w_sims[n], label="Sim $n", lw=2)
end
xlabel!("Time (101–200)")
ylabel!("Wage w")
title!("Wage paths")

plot6 = plot()
for n in 1:Nsim
    plot!(e_sims[n], label="Sim $n", lw=2)
end
xlabel!("Time (101–200)")
ylabel!("Education e")
title!("Education paths")

plot7 = plot()
for n in 1:Nsim
    plot!(y_sims[n], label="Sim $n", lw=2)
end
xlabel!("Time (101–200)")
ylabel!("Income y")
title!("Income paths")


corr_y_e = cor(y_all, e_all)

corr_w_e = cor(w_all, e_all)

#correlation is negative, that's because the cost of education grows with income (and as a result with wage shocks), and therefore we can 