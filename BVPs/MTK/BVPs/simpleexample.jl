using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables y(..)
eqs = [D(D(y(t))) ~ y(t)]

(ti, tf) = (0.0, 1.0)
cons = [y(ti) ~ 1.0, y(tf) ~ 0.0]

@named bvpsys = System(eqs, t; constraints=cons)
bvpsys = mtkcompile(bvpsys)

u0 = [0.0, 0.0] # initial guess

using InfiniteOpt
jprob = JuMPDynamicOptProblem(bvpsys, (ti, tf); dt=0.01)

#=
prob = BVProblem(bvpsys, u0, (ti, tf))
sol = solve(prob, MIRK4(), dt=0.01)
using Plots
plot(sol)
=#