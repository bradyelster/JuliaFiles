using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

# initialization of the pendulum allows for solution of initial conditions regardless of state selection

@parameters g
@variables x(t)
@variables y(t) [state_priority = 10, guess = 1.0]
@variables λ(t) [guess = 1.0]

eqns = [D(D(x)) ~ λ * x
        D(D(y)) ~ λ * y - g
        x^2 + y^2 ~ 1]

@mtkcompile pend = System(eqns, t)

# we can provide an initial condition for x and D(y)
# despite the eqns, after mtkcompile being only differential in y and D(y)

u0 = [x => 1, D(y) => 0, g => 1]
prob = ODEProblem(pend, u0, (0.0, 10.0))
sol = solve(prob, FBDF())
plot(sol; idxs = [x, y])