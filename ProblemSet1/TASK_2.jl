using LinearAlgebra, Printf

function build_matrix(α, β)
    A = zeros(5, 5)
    A[1, 1] = 1
    A[1, 2] = -1
    A[1, 4] = α - β
    A[1, 5] = β
    A[2, 2] = 1
    A[2, 3] = -1
    A[3, 3] = 1
    A[3, 4] = -1
    A[4, 4] = 1
    A[4, 5] = -1
    A[5, 5] = 1
    return A
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

# The condition number and relative residuals are higher for larger values of beta.
# x computed with the backslash operator matches the solution of the problem for these values of beta.
# Dokończyć bo chyba jeszcze trzeba wyjaśnić dlaczego tak się dzieje