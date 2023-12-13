module PointMassExamples
using Plots
using LaTeXStrings

export N_PointmassProblem,
    One_PointmassProblem,
    single_dof_damped_force,
    single_dof_damped_noforce,
    single_dof_undamped_force,
    single_dof_numerical

struct N_PointmassProblem
    ###   Constants
    k::Vector{Float64}
    m::Vector{Float64}
    γ::Vector{Float64}
    v0::Vector{Float64}
    x0::Vector{Float64}
    F::Float64
end

struct One_PointmassProblem
    ###   Constants
    k::Float64
    m::Float64
    γ::Float64
    v0::Float64
    x0::Float64
    F::Float64
end

damping_ratio(p::One_PointmassProblem) = p.γ / sqrt(2 * p.k * p.m)
nat_freq(p::One_PointmassProblem) = sqrt(p.k / p.m)

################## 
#   Analytical   #   
##################

function single_dof_damped_noforce(p::One_PointmassProblem)
    ###   2.2 Damped with no external force scenario
    ζ = damping_ratio(p)
    ω0 = nat_freq(p)

    c1 = (p.v0 + ζ * ω0 * p.x0) / (ω0 * sqrt(1 - ζ * ζ))
    c2 = x0
    a = ω0 * sqrt(1 - ζ * ζ)
    #     Solution
    x1(t) = exp(-ω0 * ζ * t) * (c1 * sin(a * t) + c2 * cos(a * t))
    #plot(x1, 0.0, 30.0, label="damped, F=0")
end


function single_dof_damped_force(p::One_PointmassProblem)
    ###   2.4 Damped with an external constant force scenario
    ζ = damping_ratio(p)
    ω0 = nat_freq(p)

    c3 = p.x0 - p.F / p.k
    c4 = (p.v0 + ω0 * ζ * c3) / (ω0 * sqrt(1 - ζ * ζ))
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
        m, c, k, f = p
        ddu[1] = -1 / m * (k * u[1] + c * du[1] + f)
    end

    prob = SecondOrderODEProblem(one_spring_damper!, p.v0, p.x0, tspan, [p.m, p.γ, p.k, p.F])
    sol = solve(prob)
end

###   Laplacian transforms
###   2.1 No damping and no external force
G(ω) = abs(1 / (-ω * ω * m + k))
plot(G, -pi, pi, xlabel="ω", ylabel="|G(ω)|", xticks=([-ω0, ω0], ["-ω0", "ω0"]))
# Think: scientific programming, asympotic behaviour

###   2.2 Damping wihtout external force
#     Real part
Hx(ω) = (-m * ω * ω + k) / ((-m * ω * ω + k)^2 + γ * γ * ω * ω)
#     Imaginary part
Hy(ω) = -(γ * ω) / ((-m * ω * ω + k)^2 + γ * γ * ω * ω)
#     Magnitude
H(ω) = sqrt(Hx(ω) * Hx(ω) + Hy(ω) * Hy(ω))

plot(H, -pi, pi, xlabel="ω", ylabel="|H(ω)|", xticks=([-ω0, ω0], ["-ω0", "ω0"]))

end