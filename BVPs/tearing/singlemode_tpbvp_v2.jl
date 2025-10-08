# Tearing Mode - Numerical Solution, 10/03/2025
# Single Helicity Mode, TwoPointBVProblem, Optimized

using BoundaryValueDiffEq, Plots, LaTeXStrings
using StaticArrays # for performance

L = 12.0
tspan = (0.0, L)
# Use @inline for small functions
@inline f(t) = tanh(t)
@inline ddf(t) = -2 * tanh(t) * sech(t)^2

k = 0.5
S = 100
ε = S^(-2 / 3)
Q_guess = 0.26706478858674904

# Pre-compute constants
const k2 = k^2
const inv_ε = 1 / ε
const inv_ε2 = 1 / ε^2

# Make function more efficient with pre-computed constants
function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    Q = p[1]

    ft = f(t)
    ddft = ddf(t)

    du[1] = ψ1
    du[2] = ϕ1
    du[3] = (Q * inv_ε + k2) * ψ - ft * inv_ε * ϕ
    du[4] = (ddft / (Q * ε) - ft * inv_ε2) * ψ + (ft^2 / (Q * ε^2) + k2) * ϕ
end

function bca!(res, u, p)
    res[1] = u[3]       # ψ'(0)=0
    res[2] = u[2]       # ϕ(0)=0
    res[3] = u[1] - 1   # ψ(0) = 1: extra constraint to fix unknown parameter Q
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

u0 = [1, 0, 0, 1]

bvp = TwoPointBVProblem(
    tearing!,
    (bca!, bcb!),
    initial_guess,
    tspan,
    [Q_guess],
    bcresid_prototype=(zeros(3), zeros(2)),
    fit_parameters=false
)

sol = solve(bvp, MIRK4(), dt=0.5)

# print the estimated value of γ which satisfies the BCs
Q = sol.prob.p[1]

scatter(sol, idxs=(0, 1), label=L"ψ(t)", continuity=:right)
scatter!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:topright, continuity=:right)
