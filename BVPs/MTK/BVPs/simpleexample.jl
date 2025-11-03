using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables y(..)
eqs = [D(D(y(t))) ~ y(t)]

@parameters deriv_guess

(ti, tf) = (0.0, 1.0)
cons = [y(ti) ~ 1.0, y(tf) ~ 0.0]

@named bvpsys = System(eqs, t; constraints=cons)
bvpsys = mtkcompile(bvpsys)

u0map = [y(t) => 1.0, D(y(t)) => deriv_guess]
parammap = [deriv_guess => 0.1]

using InfiniteOpt, Ipopt, DiffEqDevTools
jprob = JuMPDynamicOptProblem(bvpsys, [u0map; parammap], (ti, tf); dt=0.001)
jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructRadauIIA5()));

using Plots
Plots.plot(jsol.sol)

#=
prob = BVProblem(bvpsys, u0, (ti, tf))
sol = solve(prob, MIRK4(), dt=0.01)
using Plots
plot(sol)
=#