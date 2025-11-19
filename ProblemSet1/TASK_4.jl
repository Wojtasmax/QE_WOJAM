2+2
using CSV
using DataFrames
using Statistics
using Random
using Distributions
using LinearAlgebra

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

# Solve the system using the following methods:

# 1: Julia’s backslash operator

result1 = A\b
weights1 = result1[1:end-2]

# sanity check czy coś
weights1'*ones(length(weights1))
weights1'*mean_returns

# 2: 



