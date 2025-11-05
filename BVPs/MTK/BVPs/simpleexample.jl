# see my discourse post about this problem:
# https://discourse.julialang.org/t/boundary-value-problems-in-modelingtoolkit/133491

using ModelingToolkit, BoundaryValueDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables y(..)
eq = [D(D(y(t))) ~ y(t)]

# right boundary condition treated as constraint 
constr = [y(1.0) ~ 0.0]

@named bvp = System(
    eq, t;
    constraints=constr
)

bvpsys = mtkcompile(bvp)

# pass uncertain parameters/initial conditions to solver via `guesses` so it knows it can change them
prob = BVProblem(bvpsys, [y(t) => 1], (0.0, 1.0), guesses=[D(y(t)) => 0.0])
sol = solve(prob, MIRK4(), dt=0.01)
p = plot(sol)
# savefig(p, "bvp_exam.png")