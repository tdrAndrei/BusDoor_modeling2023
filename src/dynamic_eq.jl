module DynamicEq
using DifferentialEquations
using Plots
using SparseArrays
using LinearAlgebra
using Arpack
include("./BeamProblem.jl")
import .BeamProblem


export dynamic_eq!, solve_eq, plot_sol
## Constants
N = 100 #uneven number (we start counting at 1 .. index)
xp = 50
L = 1
h = L / (N - 1)
EI = 56000
load = 490
μ = 1.0
dpi = 300

A2 = BeamProblem.discretize_space(N, h)
# f = BeamProblem.point_load(N, h, xp, load)

function dynamic_eq!(ddu, du, u, p, t)
    μ, EI = p
    ddu .= (1 / μ) * (-EI * A2 * u + BeamProblem.point_load_t(N, h, xp, load, t, 0.5))
end

function solve_eq(u0, tspan)
    v0 = zeros(N)
    prob = SecondOrderODEProblem(dynamic_eq!, v0, u0, tspan, [μ, EI])
    sol = solve(prob, alg=ImplicitEuler())
    sol
end

function plot_sol(u0, tspan)
    sol = init_eq(u0, tspan)
    anim = @animate for i ∈ 1:length(sol.t)
        plot(sol.u[i][2, :], label="bending")
        ylims!(0, 0.0005)
    end every 20
    gif(anim, "dynamic_beam.gif", fps=10)
end

function eigenfrequencies()
    B = EI / μ * A2
    eigenvalues, modes = eigs(B, nev=size(B)[1])
    eigenvalues
end

end