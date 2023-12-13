module DynamicEq
using DifferentialEquations
using Plots
using SparseArrays
using LinearAlgebra
using Arpack

export BeamProblem, 
    point_load, 
    point_load_t,  
    dynamic_eq_numerical, 
    animate_sol

struct BeamProblem
    N #uneven number (we start counting at 1 .. index)
    L
    h
    A2 #What type
    EI
    μ
    xp
    u0
    F

    function BeamProblem(;N, L, μ, xp, f, EI=nothing, u0=nothing) # Constructor 
        h = L / (N - 1)
        A2 = discretize_space(N, h)
        F = f(N, h, xp)
        
        if isnothing(EI)
            EI = 10 * (L^3) / 3
        end

        if isnothing(u0)
            u0 = zeros(N)
        end

        return new(N, L, h, A2, EI, μ, xp, u0, F)
    end
end


function discretize_space(N, h)
    norm = 1 / (h * h)
    ## Compute diagonals
    v1 = -2 * norm * ones(N)
    v2 = norm * ones(N - 1)
    v3 = v2
    v1[1] = 1
    v1[end] = 1
    v2[1] = 0
    v3[end] = 0

    ## Construct the 1D forward difference matrix
    A = spdiagm(-1 => v3, 0 => v1, 1 => v2)
    ## Construct the system matrix A, [and exclude all zeros]
    A2 = dropzeros(A * A)
end

point_load(load) = (N, h, xp) -> begin
    v = zeros(N)
    v[xp] = load / h
    ## Initial conditions
    v[1] = 0
    v[end] = 0
    return (t) -> v
end

point_load_t(load, specified_t) = (N, h, xp) -> begin
    v = zeros(N)
    v[xp] = load / h
    ## Initial conditions
    v[1] = 0
    v[end] = 0
    return (t) -> begin
        if t < specified_t
            return zeros(N)
        end
        v
    end
end

function analytical_solution_static(p::BeamProblem, load)
    L, EI, h, xp = p
    # First half of the beam
    a = (xp - 1) * h
    # Last half of the beam
    b = L - a
    x = 0:h:L
    w = []
    for xi in x
        if xi <= a
            push!(w, (load * b * xi * (L * L - b * b - xi^2)) / (6 * L * EI))
        else
            push!(w, (load * b * xi * (L * L - b * b - xi^2)) / (6 * L * EI) + (load * (xi - a)^3) / (6 * EI))
        end
    end
    w
end


function dynamic_eq!(ddu, du, u, p, t)
    A2, μ, EI, F = p
    ddu .= (1 / μ) * (-EI * A2 * u + F(t))
end

function dynamic_eq_numerical(p::BeamProblem, tspan)
    v0 = zeros(p.N)
    prob = SecondOrderODEProblem(dynamic_eq!, v0, p.u0, tspan, (p.A2, p.μ, p.EI, p.F))
    sol = solve(prob, alg=ImplicitEuler())
    sol
end

function animate_sol(p::BeamProblem, tspan)
    sol = dynamic_eq_numerical(p, tspan)
    anim = @animate for i ∈ 1:length(sol.t)
        plot(sol.u[i][2, :], label="bending")
        ylims!(0, 0.0005)
    end every 20
    gif(anim, "dynamic_beam.gif", fps=10)
end

end