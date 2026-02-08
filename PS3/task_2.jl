using Distributions, Random, Statistics, Plots

# dynamic programming problem describing basil's decision problem
# basil wants the orchid to bloom before the deadline while managing stress risk

# bellman equation:
# v_t(o,s) = r(o,t) + max_a e[ v_{t+1}(o', s') ] for t < t
# v_t(o,s) = r(o,t) + φ(o,p) at the terminal period

# action set:
# a = :w (wait) or :a (apply fertilizer), only available if the orchid has not bloomed

# state variables:
# o = 0   orchid has not bloomed
# o = 1   orchid has bloomed (absorbing)
# o = -1  orchid is dead (absorbing)
# s = current stress level of the orchid

T     = 20          # finite time horizon (days until basil leaves)
Smax  = 10          # maximum admissible stress level

λ     = 0.05        # daily probability of natural bloom
δ     = 0.25        # increase in bloom probability from fertilizer
ρ     = 5.0         # poisson rate of stress increments from fertilizer

c0, c1, c2 = 0.5, 0.5, 0.1   # parameters governing anxiety costs over time
α, ψ       = 5.0, 10.0      # reward from bloom and penalty from death

p     = 100.0               # market value of the orchid
κ, θ, ω = 0.5, 2.0, 2.5     # terminal payoff coefficients

β = 1.0                     # no discounting across periods

# truncated poisson distribution for stress accumulation
# truncation keeps computations finite
Kmax = 20
pois = Poisson(ρ)
π = [pdf(pois, k) for k in 0:Kmax]

# per-period flow payoff depending on orchid state and time
function reward(O, t)
    # anxiety cost increases over time while waiting
    base_cost = c0 + c1*t + c2*t^2

    if O == 1
        # positive payoff if the orchid has bloomed
        return α
    elseif O == 0
        # negative payoff from anxiety while waiting
        return -base_cost
    else
        # additional penalty if the orchid is dead
        return -base_cost - ψ
    end
end

# terminal payoff received at time T
function terminal_payoff(O)
    if O == 1
        # basil gains utility proportional to orchid price
        return θ * p
    elseif O == 0
        # regret from never seeing the bloom
        return -κ * p
    else
        # large loss if the orchid dies
        return -ω * p
    end
end

# helper function mapping orchid states to array indices
# this keeps indexing clean and consistent
function state_index(O)
    if O == -1
        return 1
    elseif O == 0
        return 2
    else
        return 3
    end
end

# expected continuation value if basil waits
# stress does not change in this case
function EV_wait(Vnext, S)
    bloom_value    = Vnext[3, S+1]
    no_bloom_value = Vnext[2, S+1]

    return λ * bloom_value + (1 - λ) * no_bloom_value
end

# expected continuation value if basil applies fertilizer
# fertilizer increases bloom probability but raises stress stochastically
function EV_apply(Vnext, S)
    EV = 0.0

    for k in 0:Kmax
        prob = π[k+1]

        # if stress exceeds the maximum, the orchid dies
        if S + k > Smax
            EV += prob * Vnext[1, Smax+1]
        else
            bloom_value    = Vnext[3, S+k+1]
            no_bloom_value = Vnext[2, S+k+1]

            EV += prob * (
                (λ + δ) * bloom_value +
                (1 - λ - δ) * no_bloom_value
            )
        end
    end

    return EV
end

# single bellman update for a given state (t, o, s)
function bellman_update(t, O, S, V)
    oi = state_index(O)

    # absorbing states: no decisions once bloomed or dead
    if O != 0
        value = reward(O, t) + V[t+1, oi, S+1]
        return value, :None
    end

    # extract next-period value function
    Vnext = view(V, t+1, :, :)

    value_wait  = EV_wait(Vnext, S)
    value_apply = EV_apply(Vnext, S)

    # choose the action with higher expected value
    if value_apply > value_wait
        return reward(0, t) + value_apply, :A
    else
        return reward(0, t) + value_wait, :W
    end
end

# solve the model using backward induction
function solve_model()
    # value function indexed by time, orchid state, and stress
    V = zeros(T+1, 3, Smax+1)

    # policy function storing optimal actions
    policy = fill(:None, T, 3, Smax+1)

    # terminal condition
    for O in (-1, 0, 1), S in 0:Smax
        V[T+1, state_index(O), S+1] = terminal_payoff(O)
    end

    # backward induction over time and states
    for t in T:-1:1
        for O in (-1, 0, 1), S in 0:Smax
            val, act = bellman_update(t, O, S, V)
            V[t, state_index(O), S+1] = val
            policy[t, state_index(O), S+1] = act
        end
    end

    return V, policy
end

# solve the dynamic programming problem
V, policy = solve_model()

# report the optimal initial decision
println("4(a) Optimal action at t = 1, O = 0, S = 0: ", policy[1,2,1])

# plot optimal policy as a function of stress for selected times
times = [1, 5, 10, 15, 19]
plot()
for t in times
    pol = [policy[t,2,S+1] == :A ? 1 : 0 for S in 0:Smax]
    plot!(0:Smax, pol, label="t = $t", lw=2)
end
xlabel!("stress level s")
ylabel!("apply fertilizer = 1")
title!("optimal policy as a function of stress")
display(current())

# expected value at the initial state
println("4(c) Expected value V₁(0,0) = ", V[1,2,1])

# simulate paths following the optimal policy
function simulate(policy; N = 1000)
    bloom = 0
    dead  = 0
    never = 0

    fert_count = zeros(N)
    bloom_by_t = zeros(T)

    # loop over simulation runs
    for n in 1:N
        O = 0
        S = 0

        # simulate one life cycle of the orchid
        for t in 1:T

            # decision only matters if the orchid has not bloomed yet
            if O == 0
                action = policy[t, 2, S+1]

                if action == :A
                    # fertilizer is applied
                    fert_count[n] += 1

                    # stress increment from fertilizer
                    k = rand(pois)

                    if S + k > Smax
                        # orchid dies due to excessive stress
                        O = -1
                        S = Smax
                    else
                        # stress increases but remains feasible
                        S = S + k

                        # bloom may occur with higher probability
                        if rand() < λ + δ
                            O = 1
                        else
                            O = 0
                        end
                    end

                else
                    # waiting without fertilizer
                    if rand() < λ
                        O = 1
                    else
                        O = 0
                    end
                end
            end

            # once the orchid blooms, it stays bloomed forever
            if O == 1
                bloom_by_t[t:end] .+= 1
                break
            end
        end

        # record final outcome of this simulation run
        if O == 1
            bloom += 1
        elseif O == -1
            dead += 1
        else
            never += 1
        end
    end

    return bloom / N,
           never / N,
           dead / N,
           mean(fert_count),
           bloom_by_t ./ N
end

# run simulations
bloom, never, dead, avg_fert, cum_bloom = simulate(policy)

println("4(d) Fraction blooming: $bloom, never blooming: $never, dead: $dead")
println("4(e) Average number of fertilizer applications: $avg_fert")

# plot cumulative probability of bloom over time
plot(1:T, cum_bloom, lw=2)
xlabel!("time t")
ylabel!("p(oₜ = 1)")
title!("cumulative probability of bloom")
display(current())

# analyze how the value depends on the orchid price
prices = 50:10:200
Vvals = Float64[]

for ptest in prices
    global p = ptest
    Vtmp, _ = solve_model()
    push!(Vvals, Vtmp[1,2,1])
end

plot(prices, Vvals, lw=2)
xlabel!("price p")
ylabel!("expected value v₁(0,0)")
title!("expected utility as a function of orchid price")
display(current())




