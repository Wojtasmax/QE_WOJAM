using LinearAlgebra, Printf

function build_matrix(α, β)
    [
        1  -1   0     α-β   β
        0   1  -1     0     0
        0   0   1    -1     0
        0   0   0     1    -1
        0   0   0     0     1
    ]
end

function build_vector(α)
    return [α, 0, 0, 0, 1]
end

function exact_solution(α, β)
    return ones(5)
end

function solve_system(α, β)
    A = build_matrix(α, β)
    b = build_vector(α)
    x_exact = exact_solution(α, β)
    x_backslash = A \ b
    rel_res = norm(b - A * x_backslash) / norm(b)
    cond_num = cond(A)
    return (x_exact, x_backslash, rel_res, cond_num)
end

α = 0.1
β_values = 10.0 .^ (0:12)

println("Table of results for α = 0.1 and varying β:")
@printf("%12s %12s %15s %20s %20s\n", "β", "x1_exact", "x1_backslash", "condition_number", "relative_residual")
for β in β_values
    x_exact, x_backslash, rel_res, cond_num = solve_system(α, β)
    x1_exact = x_exact[1]
    x1_backslash = x_backslash[1]
    @printf("%12.1e %12.1f %15.8e %20.8e %20.8e\n", β, x1_exact, x1_backslash, cond_num, rel_res)
end

# x computed with the backslash operator matches the solution of the problem for these values of beta.
# The condition number and relative residual are higher for larger values of beta.
# High condition number expected to cause numerical instability, but doesn't in this case.
# Very small values of relative residual.
# One might expect that with large condition number values, the numerical solution would begin to deviate from the exact one, but it's possible that the optimization of the backslash operator in Julia and/or the (near-triangular) structure of the matrix A prevent this from happening.
