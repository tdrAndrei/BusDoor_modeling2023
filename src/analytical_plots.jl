using Plots

###   Constants
k     = 60
m     = 50
w0    = sqrt(k / m)
gamma = 23
zeta  = gamma / sqrt( 2 * k * m)
v0    = 0
x0    = 3
F     = 500

###   2.2 Damped with no external force scenario
c1    = (v0 + zeta * w0 * x0) / (w0 * sqrt(1 - zeta * zeta))
c2    = x0
a     = w0 * sqrt(1 - zeta * zeta)

x1(t) = exp(-w0 * zeta * t) * (c1 * sin(a * t) + c2 * cos(a * t))
plot(x1, 0.0, 30.0, label="damped, F=0")

###   2.4 Damped with an external constant force scenario
c3    = x0 - F / k
c4    = (v0 + w0 * zeta * c3) / (w0 * sqrt(1 - zeta * zeta))
b     = w0 * sqrt(1 - zeta * zeta)

x2(t) = exp(-w0 * zeta * t) * (c3 * sin(b * t) + c4 * cos(b * t)) + F / k
plot!(x2, 0.0, 30.0, label="damped, F=500")