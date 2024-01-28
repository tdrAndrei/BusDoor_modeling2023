module Beam

using DifferentialEquations
using Plots
using SparseArrays
using LinearAlgebra
using Arpack
using Distributions

export BeamProblem,
    point_load,
    point_load_t,
    dynamic_eq_numerical,
    animate_sol,
    static_eq_analytical,
    static_eq_numerical,
    gaussian_load

struct BeamProblem
    N::Int64 #uneven number (we start counting at 1 .. index)
    L::Int64
    h::Float64
    A2 #What type
    EI::Float64
    μ::Float64
    xp::Int64
    u0::Vector{Float64}
    F::Function

    function BeamProblem(; N, L, μ, xp, f, k=nothing, EI=nothing, u0=nothing, matrix_type="A2") # Constructor 
        h = L / (N - 1)
        if (matrix_type=="double_laplacian")
            A2 = double_laplacian_check(N, h, k, EI)
        elseif (matrix_type=="A2")
            A2 = discretize_space(N, h)
        end 
        
        # Default for F
        if !isnothing(f)
            F = f(N, h, xp, μ)
        else
            F = (t) -> 0
        end

        # Default for EI
        if isnothing(EI)
            EI = 10 * (L^3) / 3
        end

        #Default for u0
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

    A[1,1] = 1
    A[end,end] = 1
    ## Construct the system matrix A, [and exclude all zeros]
    A2 = dropzeros(A * A)
end

function double_laplacian(N, h, k, EI)
    norm = 1 / (h * h)
    ## Compute diagonals
    v1 = -2 * norm * ones(N)
    v2 = norm * ones(N - 1)
    v3 = norm * ones(N - 1)
    v1[1] = 1
    v1[end] = 1
    v2[1] = 0
    v3[end]=0

    v1B = -2 * norm * ones(N)
    v2B = norm * ones(N - 1)

    v2B[1] = 2 * norm 
    ## Construct the 1D forward difference matrix
    A = spdiagm(-1 => v3, 0 => v1, 1 => v2)
    B = spdiagm(-1 => v3, 0 => v1B, 1 => v2B)
    ## Construct the system matrix A, [and exclude all zeros]
    invh4 = 1/(h* h* h* h)
    A2 = A * B
    A2[end,end] =  (2 * h * h * h * k /EI + 2) * invh4
    A2[end,end-1] = -4 * invh4
    A2[end,end-2] = 2 * invh4

    A2[end-1,end-1] = 5 * invh4
    A2[end-1,end] = - 2 * invh4

    A2[1,1] = 0
    A2[1,2] = 0


    dropzeros(A2)
    
end

function double_laplacian_check(N, h, k, EI)
    norm = 1 / (h * h * h * h)
    ## Compute diagonals
    v1 = 1 * norm * ones(N-2)
    v2 = -4 * norm * ones(N-1)
    v3 = 6 * norm * ones(N)
    v4 = -4 * norm * ones(N-1)
    v5 = 1 * norm * ones(N-2)

    v3[2] = 5 * norm
    v3[1] = 1 * norm
    v3[end-1] = 5 * norm
    v4[1] = 0 * norm
    v5[1] = 0 * norm


    ## Construct the 1D forward difference matrix
    C = spdiagm(-2 => v1, -1 => v2, 0 => v3, 1 => v4, 2 => v5)
    ## Construct the system matrix A, [and exclude all zeros]
    C[end,end] = (2 * h * h * h * k /EI + 2) * norm
    C[end,end-1] = -4 * norm
    C[end,end-2] = 2 * norm

    C[end-1,end] = -2 * norm
    C[end-1,end-1] = 5 * norm

    dropzeros(C)
end

################################
#            Forces            #
################################

point_load(load) = (N, h, xp, μ) -> begin
    v = zeros(N)
    v[xp] = load / h
    return (t) -> v
end

point_load_t(load, specified_t) = (N, h, xp, μ) -> begin
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

gaussian_load(load) = (N, h, xp, μ) -> begin
    if (xp > N)
        print("index too large")
        return
    end

    f = zeros(N)

    n = 30
    c = xp*h
    σ = 0.05

    last = min(xp+n, N)  # Ensure last does not exceed N
    points = xp+n > N ? (2 * n - (xp+n - N) + 1) : (2 * n + 1)

    d = Normal(c, σ)

    # Calculate quantiles
    lo = max((xp - n) * h, -5σ + c)
    hi = min((last) * h, 5σ + c)

    # Generate x values
    x = range(lo, hi; length = points)

    # Calculate PDF values for the x values
    pdf_values = pdf.(d, x)

    area = h * sum(pdf_values)
    pdf_values .= pdf_values * load/(area)

    f[xp-n:last] = pdf_values
    return (t) -> f
end


################################
#       Static Equations       #
################################

function static_eq_analytical(p::BeamProblem, load)
    (; L, EI, h, xp) = p
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

function static_eq_numerical(p::BeamProblem)
    (; A2, EI, F) = p
    qx = F(0) #F is constant in time, so get the value at any point
    return A2 \ (qx / EI)
end

#################################
#       Dynamic Equations       #
#################################

function dynamic_eq!(ddu, du, u, p, t)
    A2, μ, EI, F = p
    ddu .= (1 / μ) * (-EI * A2 * u + F(t))
    #dddu .= 0
end

function dynamic_eq_numerical(p::BeamProblem, tspan)
    v0 = zeros(p.N)
    prob = SecondOrderODEProblem(dynamic_eq!, v0, p.u0, tspan, (p.A2, p.μ, p.EI, p.F))
    sol = solve(prob, alg=ImplicitEuler())
    sol
end

function animate_sol(sol, h, L, ymax=0.0005)
    anim = @animate for i ∈ 1:length(sol.t)
        plot((0:h:L),sol.u[i][2, :], xlabel="Length[m]",ylabel="Displacement[m]",label="bending t=$(sol.t[i])")
        ylims!(0, ymax)
    end every 20
    gif(anim, "dynamic_beam.gif", fps=10)
end

end
# Beam.discretize_space(5, 1)
# Beam.double_laplacian(6, 1, 3, 1)
# Beam.double_laplacian_check(7, 1, 3, 1)