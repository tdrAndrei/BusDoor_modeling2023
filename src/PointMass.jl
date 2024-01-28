module PointMass

using Plots
using LaTeXStrings
using DifferentialEquations

export N_PointmassProblem,
    One_PointmassProblem,
    single_dof_damped_force,
    single_dof_damped_noforce,
    single_dof_undamped_force,
    single_dof_numerical,
    n_dof_numerical

################## 
#      Types     #   
##################

struct N_PointmassProblem
    ###   Constants
    m::Vector{Float64}
    γ::Vector{Float64}
    k::Vector{Float64}
    v0::Vector{Float64}
    x0::Vector{Float64}
    F::Vector{Float64}

    ###   Computed
    M::Matrix{Float64}
    C::Matrix{Float64}
    K::Matrix{Float64}
    function N_PointmassProblem(; m::Vector{Float64}, γ::Vector{Float64}, k::Vector{Float64}, v0::Vector{Float64}, x0::Vector{Float64}, F::Vector{Float64})
        n = size(m, 1)
        M = zeros((n, n))
        K = zeros((n, n))
        C = zeros((n, n))

        for i in 1:n
            M[i, i] = m[i]
            K[i, i] = k[i] + k[i+1]
            C[i, i] = γ[i] + γ[i+1]

            if i > 1
                K[i, i-1] = -k[i]
                C[i, i-1] = -γ[i]
            end

            if i < n
                K[i, i+1] = -k[i+1]
                C[i, i+1] = -γ[i+1]
            end
        end

        new(m, γ, k, v0, x0, F, M, C, K)
    end

end


Base.@kwdef struct One_PointmassProblem
    ###   Constants
    m::Float64
    γ::Float64
    k::Float64
    v0::Float64
    x0::Float64
    F::Float64
end

damping_ratio(p::One_PointmassProblem) = p.γ / (2 * sqrt(p.k * p.m))
nat_freq(p::One_PointmassProblem) = sqrt(p.k / p.m)

################## 
#   Analytical   #   
##################

function single_dof_damped_noforce(p::One_PointmassProblem)
    ###   2.2 Damped with no external force scenario
    ζ = damping_ratio(p)
    ω0 = nat_freq(p)

    c1 = (p.v0 + ζ * ω0 * p.x0) / (ω0 * sqrt(1 - ζ * ζ))
    c2 = p.x0
    a = ω0 * sqrt(1 - ζ * ζ)
    #     Solution
    x1(t) = exp(-ω0 * ζ * t) * (c1 * sin(a * t) + c2 * cos(a * t))
    #plot(x1, 0.0, 30.0, label="damped, F=0")
end


function single_dof_damped_force(p::One_PointmassProblem)
    ###   2.4 Damped with an external constant force scenario
    ζ = damping_ratio(p)
    ω0 = nat_freq(p)

    c4 = p.x0 - p.F / p.k
    c3 = (p.v0 + ω0 * ζ * c4) / (ω0 * sqrt(1 - ζ * ζ))
    b = ω0 * sqrt(1 - ζ * ζ)
    #     Solution
    x2(t) = exp(-ω0 * ζ * t) * (c3 * sin(b * t) + c4 * cos(b * t)) + p.F / p.k
    #plot!(x2, 0.0, 30.0, label="damped, F=500")
end


function single_dof_undamped_force(p::One_PointmassProblem)
    ###   2.3 Not damped with an external constant force scenario
    ω0 = nat_freq(p)

    c5 = p.x0 - p.F / p.k
    c6 = p.v0 / ω0
    #     Solution
    x3(t) = c5 * cos(ω0 * t) + c6 * sin(ω0 * t)
    #plot!(x3, 0.0, 30.0, label="not damped, F=500")
end

################## 
#   Numerical    #   
##################

function single_dof_numerical(p::One_PointmassProblem, tspan)
    function one_spring_damper!(ddu, du, u, p, t)
        m, γ, k, f = p
        ddu[1] = -1 / m * (k * u[1] + γ * du[1] - f)
    end

    prob = SecondOrderODEProblem(one_spring_damper!, [p.v0], [p.x0], tspan, [p.m, p.γ, p.k, p.F])
    sol = solve(prob)
end

function n_dof_numerical(p::N_PointmassProblem, tspan)
    function N_spring_damper!(ddu, du, u, p, t)
        M, C, K, F = p
        ddu .= M \ (-C * du - K * u + F)
    end

    prob = SecondOrderODEProblem(N_spring_damper!, p.v0, p.x0, tspan, [p.M, p.C, p.K, p.F])
    sol = solve(prob)
end

########################## 
#   Laplace transforms   #   
##########################

function single_dof_undamped_noforce_laplace(p::One_PointmassProblem)
    ###   2.1 No damping and no external force
    ω0 = nat_freq(p)
    G(ω) = abs(1 / (-ω * ω * p.m + p.k))
    plot(G, -pi, pi, xlabel="ω", ylabel="|G(ω)|", xticks=([-ω0, ω0], ["-ω0", "ω0"]))
    # Think: scientific programming, asympotic behaviour
end

function single_dof_damped_noforce_laplace(p::One_PointmassProblem)
    ###   2.2 Damping wihtout external force
    ω0 = nat_freq(p)
    #     Magnitude
    H(ω) = abs(1 / (m * (im * ω)^2 + γ * (im * ω) + k))

    plot(H, -pi, pi, xlabel="ω", ylabel="|H(ω)|", xticks=([-ω0, ω0], ["-ω0", "ω0"]))
end

end