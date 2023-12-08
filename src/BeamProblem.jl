module BeamProblem
using SparseArrays

export discretize_space, point_load, analytical_solution_static

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

function point_load(N, h, xp, load)
    v = zeros(N)
    v[xp] = load / h
    ## Initial conditions
    v[1] = 0
    v[end] = 0
    v
end

function point_load_t(N, h, xp, load, t, specified_t)
    v = zeros(N)
    if (t > specified_t)
        return v
    end

    v[xp] = load / h
    ## Initial conditions
    v[1] = 0
    v[end] = 0
    v
end

function analytical_solution_static(N, h, L, x, xp, load, EI)
    # First half of the beam
    a = (xp - 1) * h
    # Last half of the beam
    b = L - a
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

end