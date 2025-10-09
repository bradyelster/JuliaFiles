# Let's solve the simple Van der Pol oscillator as an IVP
# we will use the solution to solve again as a TwoPointBVP
# and we will use an implicit solver and give it the Jacobian (analytic)

using OrdinaryDiffEq, Plots

ε = 0.1
tspan = (0, 20)

function vdp!(du, u, p, t)
    x, x1 = u
    du[1] = x1
    du[2] = -x + ε * (1 - x^2) * x1
end

u0 = [1, 0]

prob = ODEProblem(vdp!, u0, tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

plot(sol)

sol.u[end]