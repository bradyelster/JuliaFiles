# This example is from Strogatz's online course on Asymptotic Methods
# found here: https://www.youtube.com/watch?v=UiO0Y_eFVQY&list=PL5EH0ZJ7V0jV7kMYvPcZ7F9oaf_YAlfbI&index=13

using BoundaryValueDiffEq, Plots, BenchmarkTools

#=
# the following doesn't do anything?
using ProgressLogging, TerminalLoggers
import Logging: global_logger
import TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
=#

function ode!(du, u, p, t)
    du[1] = u[2]
    du[2] = -(1 / p[1] + 1) * u[2] - (1 / p[1]) * u[1]
end

# Boundary conditions
function bca!(res, ua, p)
    res[1] = ua[1] - 0
end

function bcb!(res, ub, p)
    res[1] = ub[1] - 1
end

# Time span
tspan = (0.0, 1.0)

# Starting parameters
ε = 0.1
sol = [0.0, 1.0]

# Continuation loop
for i in 1:2
    global eps = ε / 10
    prob = TwoPointBVProblem(ode!, (bca!, bcb!), sol, tspan, [eps],
        bcresid_prototype=(zeros(1), zeros(1)))
    global sol = solve(prob, MIRK4(), dt=0.05, reltol=1e-9, abstol=1e-10, adaptive=true, verbose=true, progress=true)
end

plot(sol, idxs=(0, 1), label="ε = $eps", lw=2, xlabel="x", ylabel="y(x)", legend=:best, xlims=(0, 50 * eps))
vline!([ε], lw=2, linestyle=:dash, label="δ")