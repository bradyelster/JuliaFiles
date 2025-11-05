using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D
using BoundaryValueDiffEq

@variables x(..) = 1 v(..) = 0
@variables u(..) = 1 [bounds = (-1.0, 1.0), input = true]
constr = [v(1.0) ~ 0.0] # constraint
cost = [-x(t)] # cost to maximize x(t) (minimize -x(t))
@named block = System(
    [D(x(t)) ~ v(t), D(v(t)) ~ u(t)], t;
    constraints=constr
)

blocksys = mtkcompile(block; inputs=[u(t)])

#u0 = [x => 1, v => 0, u => 1]
prob = BVProblem(blocksys, [], (0.0, 1.0))
sol = solve(prob, MIRK4(), dt=0.5)
plot(sol)