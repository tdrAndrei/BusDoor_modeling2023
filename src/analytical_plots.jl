using Plots
using LaTeXStrings

###   Constants
k     = 60
m     = 50
ω0    = sqrt(k / m)
γ     = 23
ζ     = γ / sqrt( 2 * k * m)
v0    = 0
x0    = 3
F     = 500

###   2.2 Damped with no external force scenario
c1    = (v0 + ζ * ω0 * x0) / (ω0 * sqrt(1 - ζ * ζ))
c2    = x0
a     = ω0 * sqrt(1 - ζ * ζ)
#     Solution
x1(t) = exp(-ω0 * ζ * t) * (c1 * sin(a * t) + c2 * cos(a * t))
plot(x1, 0.0, 30.0, label="damped, F=0")


###   2.4 Damped with an external constant force scenario
c3    = x0 - F / k
c4    = (v0 + ω0 * ζ * c3) / (ω0 * sqrt(1 - ζ * ζ))
b     = ω0 * sqrt(1 - ζ * ζ)
#     Solution
x2(t) = exp(-ω0 * ζ * t) * (c3 * sin(b * t) + c4 * cos(b * t)) + F / k
plot!(x2, 0.0, 30.0, label="damped, F=500")

###   2.3 Not damped with an external constant force scenario
c5    = x0 - F / k
c6    = v0 / ω0
#     Solution
x3(t) = c5 * cos(ω0 * t) + c6 * sin(ω0 * t)
plot!(x3, 0.0, 30.0, label="not damped, F=500")