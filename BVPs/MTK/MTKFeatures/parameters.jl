using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t)
@parameters f(::Real)::Real q

@mtkcompile sys = System([D(x) ~ f(t) + q * x], t)

prob = ODEProblem(sys, [x => 1.0, f => sin, q => -1.0], (0.0, 5.0))
sol = solve(prob, Tsit5())

theme(:juno)
plot(sol)

prob.ps[f] => sin
sol = solve(prob, Tsit5())
plot(sol)