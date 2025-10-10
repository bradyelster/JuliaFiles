# Tearing Mode - Numerical Solution, 10/03/2025
# Single Helicity Mode, TwoPointBVProblem, Optimized

using BoundaryValueDiffEq, Plots, LaTeXStrings
using StaticArrays # for performance

L = 15.0
tspan = (0.0, L)
# Use @inline for small functions
@inline f(t) = tanh(t)

# m = 2
k = 0.5
# k0 = 0.25
S = 1000
# ε = S^(-2 / 3)
γ_guess = 0.001325864736286539 # 0.025988534151142495 for non-resonant mode @ S=100, k = 0.5 
# 0.26706478858674904*S^(-1/3) for resonant mode with different normalization

function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    γ = p[1]

    du[1] = ψ1
    du[2] = ϕ1
    du[3] = (-S * k * f(t)) * ϕ + (S * γ - k^2) * ψ
    du[4] = ((k^2 * S * f(t)^2) / γ) * ϕ - (k * f(t) / γ) * (-2 * k^2 + S*γ + 2 * sech(t)^2) * ψ
end

function bca!(res, u, p)
    res[1] = u[3]       # ψ'(0)=0
    res[2] = u[2]       # ϕ(0)=0
    res[3] = u[1] - 1   # ψ(0) = 1: extra constraint to fix unknown parameter
end

function bcb!(res, u, p)
    res[1] = u[1]       # ψ(L)=0
    res[2] = u[2]       # ϕ(L)=0
end

# Optimize initial guess vector
function initial_guess(p, t)
    exp_t2 = exp(-t^2)
    t_exp = t * exp_t2
    SVector(exp_t2, t_exp, -2 * t_exp, (1 - 2 * t^2) * exp_t2)
end

bvp = TwoPointBVProblem(
    tearing!,
    (bca!, bcb!),
    initial_guess,
    tspan,
    [γ_guess],
    bcresid_prototype=(zeros(3), zeros(2)),
    fit_parameters=true
)

sol = solve(bvp, MIRK4(), dt=0.01, saveat=0.5, verbose=true)

# print the estimated value of γ which satisfies the BCs
sol.prob.p[1]

plot(sol, idxs=(0, 1), label=L"ψ(t)", continuity=:right, lw=2)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:topright, continuity=:right, lw=2)