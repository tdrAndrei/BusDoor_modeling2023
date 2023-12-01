using SparseArrays
using Plots
include("./BeamProblem.jl")
import .BeamProblem

## Constants
N = 100; #uneven number (we start counting at 1 .. index)
xp = 50;
L = 1;
h = L / (N - 1);
EI = 56000;
load = 490;
dpi = 300

## 
A2 = BeamProblem.discretize_space(N, h);
qx = BeamProblem.point_load(N, h, xp, load);

## Solutions to plot
x = 0:h:L
yA = BeamProblem.analytical_solution_static(N, h, L, x, xp, load, EI)
yN = A2 \ (qx / EI)                 ## Iterative solver? Equivalent of spsolve in Python?

## Plotting
plot(x, yN, label="Numerical Solution", xlabel="Length [m]", ylabel="Deflection [m]", show=true, dpi=dpi)
plot!(x, yA, label="Analytical Solution", linestyle=:dash, color=:red, show=true, dpi=dpi)