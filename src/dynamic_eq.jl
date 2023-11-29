using DifferentialEquations
using Plots
using SparseArrays
include("./BeamProblem.jl")
import .BeamProblem

## Constants
N = 100; #uneven number (we start counting at 1 .. index)
xp = 50;
L = 1;
h = L / (N - 1);
EI = 1.0;
load = 490;
μ = 1.0
dpi = 300

A2 = BeamProblem.discretize_space(N, h)

function dynamic_eq!(ddu, du, u, p, t)
    μ, EI = p
    ddu .= (1 / μ) * (-EI * A2 * u + BeamProblem.point_load(N, h, xp, load))
end

function init_eq(u0, tspan)
    v0 = zeros(N)
    prob = SecondOrderODEProblem(dynamic_eq!, v0, u0, tspan, [μ, EI])
    sol = solve(prob)
    sol
end

u0 = zeros(N)
tspan = (0.0, 10.0)
sol = init_eq(u0, tspan)
anim = @animate for i ∈ 1:length(sol.t)
    plot(sol.u[i], label="bending")
end every 5

gif(anim, "dynamic_beam.gif", fps=10)