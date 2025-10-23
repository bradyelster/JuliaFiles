using BoundaryValueDiffEq
using Plots
function f!(du, u, p, t)
    e = 2.7
    du[1] = (1 + u[2] - sin(t)) * u[4] + cos(t)
    du[2] = cos(t)
    du[3] = u[4]
    du[4] = (u[1] - sin(t)) * (u[4] - e^t)
end

function bc!(res, u, p, t)
    res[1] = u(0.0)[1]
    res[2] = u(0.0)[3] - 1
    res[3] = u(1.0)[2] - sin(1.0)
end

u0 = [0.0, 0.0, 0.0, 0.0]
tspan = (0.0, 1.0)

fun = BVPFunction(f!, bc!, mass_matrix=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0])

prob = BVProblem(fun, u0, tspan)

sol = solve(prob, RadauIIa7(), abstol=1e-10, reltol=1e-9, dt=0.05)

plot(sol)