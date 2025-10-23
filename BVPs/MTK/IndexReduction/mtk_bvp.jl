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

#=
sol = solve(prob,
    Ascher4(; zeta=[0.0, 1.0, 0.0], jac_alg=BVPJacobianAlgorithm(AutoForwardDiff()));
    dt=0.01)
plot(sol)
=#

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel dae_example begin
    @constants begin
        ε = 1
    end
    @variables begin
        x1(t) = 0
        x2(t) = 0
        x3(t) = 0
        y(t) = 0
        p1(t)
        p2(t)
    end
    @equations begin
        p1 ~ sin(t)
        p2 ~ sin(t)
        D(x1) ~ (ε + x2 - p2) * y + D(p1)
        D(x2) ~ D(p2)
        D(x3) ~ y
        0 ~ (x1 - p1) * (y - exp(t))
    end
end

@mtkcompile sys = dae_example()
@show equations(sys)