# Let's solve the simple Van der Pol oscillator as TwoPointBVP
# and we will use an implicit solver and give it the Jacobian (analytic)

using BoundaryValueDiffEq, Plots

ε = 0.1
tspan = (0, 20)

function vdp!(du, u, p, t)
    x, x1 = u
    du[1] = x1
    du[2] = -x + ε * (1 - x^2) * x1
end

function vdp_jac!(J, u, p, t)
    ε = p[1]
    x, x1 = u
    J[1, 1] = 0
    J[1, 2] = 1
    J[2, 1] = -2 * x * x1 * ε - 1
    J[2, 2] = ε * (1 - x^2)
end

function bc!(res, u, p, t)
    res[1] = u[1][1] - 1.0  # Left boundary condition: x(0) = 1
    res[2] = u[end][1] - 0.649067888450811  # Right boundary condition: x(20) = 0.649...
end

u0 = [1.0, 0.0]

bvp_f = BVPFunction(vdp!, bc!; jac=vdp_jac!, bcresid_prototype=zeros(2))
prob = BVProblem(bvp_f, u0, tspan)
@time sol = solve(prob, MIRK4(), dt=0.05, reltol=1e-8, abstol=1e-8)

plot(sol)

sol.u[end]