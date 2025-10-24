# Example of solving a BVP, DAE using ModelingToolkit 
# difference between numeric and symbolic solutions is zero!

using BoundaryValueDiffEq, Plots

p1(t) = sin(t)
dp1(t) = cos(t)
p2(t) = sin(t)
dp2(t) = cos(t)

function f!(du, u, p, t)
    x1, x2, x3, y = u
    ε = p[1]
    du[1] = (ε + x2 - p2(t)) * y + dp1(t)
    du[2] = dp2(t)
    du[3] = y
    du[4] = (u[1] - p1(t)) * (y - exp(t))
end

function bc!(res, u, p, t)
    res[1] = u[1] - p1(0)
    res[2] = u[2] - p2(1.0)
    res[3] = u[3] - 1.0
end

u0 = [0.0, 0.0, 0.0, 0.0]

tspan = (0.0, 1.0)

fun = BVPFunction(f!, bc!, mass_matrix=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0])

prob = BVProblem(fun, u0, tspan, [1])

sol = solve(prob,
    Ascher4(; zeta=[0.0, 1.0, 0.0], jac_alg=BVPJacobianAlgorithm(AutoForwardDiff()));
    dt=0.01)
plot(sol, labels=["x₁" "x₂" "x₃" "y"])

using Printf

@printf "x1(0) = %.6f (expected %.6f)\n" sol(0)[1] sin(0)
@printf "x2(1) = %.6f (expected %.6f)\n" sol(1)[2] sin(1)
@printf "x3(0) = %.6f (expected %.6f)\n" sol(0)[3] 1.0

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel dae_example begin
    @constants begin
        ε = 1
    end
    @variables begin
        # initialize with only BCs we know at t=0
        # Recall the BCs:
        # x1(0) = p1(0)
        # x2(1) = p2(1)
        # x3(0) = 1
        p1(t)
        p2(t)
        x1(t) = p1
        x2(t) = 0
        x3(t) = 1
        y(t) = 0
    end
    @equations begin
        # DAE equations
        p1 ~ sin(t)
        p2 ~ sin(t)
        D(x1) ~ (ε + x2 - p2) * y + D(p1)
        D(x2) ~ D(p2)
        D(x3) ~ y
        0 ~ (x1 - p1) * (y - exp(t))
    end
end

@mtkcompile sys = dae_example()
mtkprob = BVProblem(sys, [], tspan)
mtksol = solve(mtkprob, MIRK4(), dt=0.01)
plot(mtksol)
@show unknowns(sys)

@printf "x1(0) = %.6f (expected %.6f)\n" mtksol(0)[4] sin(0)
@printf "x2(1) = %.6f (expected %.6f)\n" mtksol(1)[3] sin(1)
@printf "x3(0) = %.6f (expected %.6f)\n" mtksol(0)[2] 1.0