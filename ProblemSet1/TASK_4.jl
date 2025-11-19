2+2
using CSV
using DataFrames
using Statistics
using Random
using Distributions
using LinearAlgebra
using IterativeSolvers
#1. load data into a df
cd(@__DIR__)
df = CSV.read("asset_returns.csv", DataFrame)
returns=Matrix(df)

# Compute the sample mean vector µ and sample covariance matrix Σ of asset returns

mean_returns = vec(mean(returns, dims=1))
cov_matrix = cov(returns)

# Form the (n + 2) × (n + 2) matrix system Ax = x above with target return µ¯ = 0.10 
μ=.1

A = vcat(
    hcat(cov_matrix, mean_returns, ones(size(cov_matrix, 1), 1)),
    hcat(mean_returns', 0.0, 0.0),
    hcat(ones(1, size(cov_matrix, 2)), 0.0, 0.0)
    )

b=vcat(zeros(size(cov_matrix, 1)), μ, 1.0)

# Report the condition number of the A matrix
cond_A = cond(A)

# custom function for checking if conditions are satisfied
function sanity_check(weights)
    sum_weights = sum(weights)
    portfolio_return = dot(weights, mean_returns)
    return sum_weights, portfolio_return
end 

# Solve the system using the following methods:

# 1: Julia’s backslash operator

result1 = A\b
weights1 = result1[1:end-2]

# sanity check czy coś
sanity_check(weights1)

# 2 Jacobi or Gauss-Seidel: 

gramA = A'*A
#check for zeros on diagonal
any(diag(gramA) .== 0)
#false so there aint no zeros, great
cond_gramA = cond(gramA)
# TODO check diagonal dominance (by columns)
function check_diag_dominance(Mat::Matrix)
    # 1- diag dominant 0- guess what
    for i in axes(Mat, 1)
        if 2*Mat[i, i]<sum(Mat[:,i])
            return 0
        end
    end
    return 1
end

check_diag_dominance(gramA)

# It seems that the A^T A matrix is not diagonally dominant (it doesnt look so anyway)
# So there is no need to write this method indeed the custom one returns nonsense
# TODO write custom, for now I use the built in one
result2=jacobi(gramA, A'*b)
weights2=result2[1:end-2]
sanity_check(weights2)

# 3 Conjugate Gradient method (CG) 
issymmetric(gramA)
isposdef(gramA)
result3=cg(gramA, A'*b)
weights3=result3[1:end-2]
sanity_check(weights3)

# 4 GMRES (the generalized minimal residual method)
result4=gmres(A, b)
weights4=result4[1:end-2]
sanity_check(weights4)

# 5 sleep
