#!/usr/bin/env julia

# Stationary-variance Lyapunov benchmark.
#
# Solves
#
#   A * Sigma + Sigma * A' = V
#
# for the stationary variance Sigma.  Three paths are timed:
#
#   1. dense Julia lyap(A, -V), which computes Schur internally;
#   2. planned dense Schur: precompute A = Z*T*Z' once, then time
#      Z'VZ, triangular/quasi-triangular LAPACK trsyl!, and ZYZ';
#   3. SMBP with rotation: rotate V into the 2x2 block basis, solve each
#      block-pair Lyapunov system using the equal-diagonal formulas used by
#      BlockDiagonalLyapunovSolver, then rotate Sigma back.

using LinearAlgebra
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "schur_vs_smbp_exp_adjoint.jl"))
include(joinpath(@__DIR__, "smbp_lyapunov_plan.jl"))

const LYAP_DEFAULT_K_VALUES = [4, 8, 16, 32, 64]

struct LyapunovInstance
    K::Int
    R::Matrix{Float64}
    D::Matrix{Float64}
    diag::Vector{Float64}
    upper::Vector{Float64}
    lower::Vector{Float64}
    A::Matrix{Float64}
    V_block::Matrix{Float64}
    V_dense::Matrix{Float64}
end

mutable struct DenseSchurLyapunovPlan
    K::Int
    Z::Matrix{Float64}
    T::Matrix{Float64}
    W::Matrix{Float64}
    tmp::Matrix{Float64}
    out::Matrix{Float64}
end

function DenseSchurLyapunovPlan(A::Matrix{Float64})
    factor = schur(A)
    K = size(A, 1)
    return DenseSchurLyapunovPlan(
        K,
        Matrix(factor.Z),
        Matrix(factor.T),
        zeros(K, K),
        zeros(K, K),
        zeros(K, K),
    )
end

mutable struct SMBPLyapunovWorkspace
    V_block::Matrix{Float64}
    Sigma_block::Matrix{Float64}
    tmp::Matrix{Float64}
    out::Matrix{Float64}
end

SMBPLyapunovWorkspace(K::Int) = SMBPLyapunovWorkspace(
    zeros(K, K),
    zeros(K, K),
    zeros(K, K),
    zeros(K, K),
)

function make_lyapunov_instance(rng::AbstractRNG, K::Int)
    smbp = make_instance(rng, K)
    G = randn(rng, K, K)
    V_block = Symmetric((G' * G) / K + 0.50I) |> Matrix
    V_dense = smbp.R * V_block * smbp.R'
    return LyapunovInstance(
        K,
        smbp.R,
        smbp.D,
        smbp.diag,
        smbp.upper,
        smbp.lower,
        smbp.A,
        V_block,
        V_dense,
    )
end

function symmetrize_matrix!(A::Matrix{Float64})
    K = size(A, 1)
    @inbounds for j in 1:K
        for i in (j + 1):K
            v = 0.5 * (A[i, j] + A[j, i])
            A[i, j] = v
            A[j, i] = v
        end
    end
    return A
end

function dense_lyapunov(A::Matrix{Float64}, V::Matrix{Float64})
    X = lyap(A, -V)
    return symmetrize_matrix!(X)
end

function dense_schur_lyapunov!(plan::DenseSchurLyapunovPlan, V::Matrix{Float64})
    mul!(plan.tmp, plan.Z', V)
    mul!(plan.W, plan.tmp, plan.Z)
    _, scale = LAPACK.trsyl!('N', 'T', plan.T, plan.T, plan.W, 1)
    rmul!(plan.W, inv(scale))
    mul!(plan.tmp, plan.Z, plan.W)
    mul!(plan.out, plan.tmp, plan.Z')
    return symmetrize_matrix!(plan.out)
end

function smbp_lyapunov_with_rotation!(work::SMBPLyapunovWorkspace,
                                      plan::SMBPLyapunovEqualDiagonalPlan,
                                      inst::LyapunovInstance,
                                      V_dense::Matrix{Float64})
    mul!(work.tmp, inst.R', V_dense)
    mul!(work.V_block, work.tmp, inst.R)
    symmetrize_matrix!(work.V_block)
    apply_smbp_lyapunov_symmetric!(work.Sigma_block, plan, work.V_block)
    mul!(work.tmp, inst.R, work.Sigma_block)
    mul!(work.out, work.tmp, inst.R')
    return symmetrize_matrix!(work.out)
end

function run_lyapunov_case(inst::LyapunovInstance; warmup::Int, reps::Int, batch::Int)
    K = inst.K
    schur_plan = DenseSchurLyapunovPlan(inst.A)
    smbp_plan = SMBPLyapunovEqualDiagonalPlan(inst.diag, inst.upper, inst.lower)
    smbp_work = SMBPLyapunovWorkspace(K)

    Sigma_dense = dense_lyapunov(inst.A, inst.V_dense)
    Sigma_schur = dense_schur_lyapunov!(schur_plan, inst.V_dense)
    Sigma_smbp = smbp_lyapunov_with_rotation!(smbp_work, smbp_plan, inst, inst.V_dense)

    schur_relerr = relative_error(Sigma_dense, Sigma_schur)
    smbp_relerr = relative_error(Sigma_dense, Sigma_smbp)
    residual_relerr = relative_error(inst.A * Sigma_smbp + Sigma_smbp * inst.A', inst.V_dense)

    dense_stats = median_time_ms(
        () -> begin
            Sigma = dense_lyapunov(inst.A, inst.V_dense)
            return Sigma[1, 1]
        end;
        warmup = warmup, reps = reps, batch = batch)

    schur_stats = median_time_ms(
        () -> begin
            Sigma = dense_schur_lyapunov!(schur_plan, inst.V_dense)
            return Sigma[1, 1]
        end;
        warmup = warmup, reps = reps, batch = batch)

    smbp_stats = median_time_ms(
        () -> begin
            Sigma = smbp_lyapunov_with_rotation!(smbp_work, smbp_plan, inst, inst.V_dense)
            return Sigma[1, 1]
        end;
        warmup = warmup, reps = reps, batch = batch)

    return (
        K = K,
        blocks = div(K, 2),
        reps = reps,
        batch = batch,
        schur_relerr = schur_relerr,
        smbp_relerr = smbp_relerr,
        residual_relerr = residual_relerr,
        dense_lyap_ms = dense_stats[1],
        schur_known_lyap_ms = schur_stats[1],
        smbp_rotation_lyap_ms = smbp_stats[1],
    )
end

function lyapunov_header()
    return "K,blocks,reps,batch,schur_relerr,smbp_relerr,residual_relerr," *
           "dense_lyap_ms,schur_known_lyap_ms,smbp_rotation_lyap_ms," *
           "dense_over_schur_known,dense_over_smbp_rotation," *
           "schur_known_over_smbp_rotation"
end

function write_lyapunov_row(io, row)
    values = Any[
        row.K, row.blocks, row.reps, row.batch,
        row.schur_relerr, row.smbp_relerr, row.residual_relerr,
        row.dense_lyap_ms, row.schur_known_lyap_ms, row.smbp_rotation_lyap_ms,
        row.dense_lyap_ms / row.schur_known_lyap_ms,
        row.dense_lyap_ms / row.smbp_rotation_lyap_ms,
        row.schur_known_lyap_ms / row.smbp_rotation_lyap_ms,
    ]
    println(io, join(fmt_csv.(values), ","))
end

function main_lyapunov()
    k_values = parse_int_list(get(ENV, "K_VALUES", nothing), LYAP_DEFAULT_K_VALUES)
    warmup = env_int("WARMUP", 3)
    base_reps = env_int("REPS", 10)
    batch = env_int("INNER_REPS", 50)
    seed = env_int("SEED", 17)
    out_path = get(ENV, "OUT", joinpath(@__DIR__, "lyapunov_stationary_variance_results.csv"))

    rng = MersenneTwister(seed)
    rows = []
    for K in k_values
        inst = make_lyapunov_instance(rng, K)
        reps = reps_for_K(K, base_reps)
        @info "running Lyapunov case" K reps batch
        row = run_lyapunov_case(inst; warmup = warmup, reps = reps, batch = batch)
        push!(rows, row)
        @printf("K=%d err_schur=%.3e err_smbp=%.3e dense=%.4fms schur=%.4fms smbp=%.4fms\n",
                row.K, row.schur_relerr, row.smbp_relerr,
                row.dense_lyap_ms, row.schur_known_lyap_ms, row.smbp_rotation_lyap_ms)
    end

    open(out_path, "w") do io
        println(io, lyapunov_header())
        for row in rows
            write_lyapunov_row(io, row)
        end
    end
    println("wrote ", out_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_lyapunov()
end
