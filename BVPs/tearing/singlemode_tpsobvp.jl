# Tearing Mode - Numerical Solution, 10/01/2025
# Single Helicity Mode, TwoPointSecondOrderBVProblem

using BoundaryValueDiffEq, Plots, LaTeXStrings

L = 12.0
tspan = (0.0, L)
f(t) = tanh(t)
ddf(t) = -2 * tanh(t) * sech(t)^2

k = 0.5
S = 100
ε = S^(-2/3)
Q_guess = 0.26706479450097753 # Q -> 1 is the NCF limit where γ = S^(-1/3), CF limit is Q << 1

function tearing!(ddu, du, u, p, t)
    ψ, ϕ = u
    Q = p[1]
    ddu[1] = (Q/ε + k^2)*ψ - (f(t)/ε)*ϕ
    ddu[2] = (ddf(t)/(Q*ε) - f(t)/ε^2)*ψ + (f(t)^2/(Q*ε^2) + k^2)*ϕ
end

function bca!(res, du, u, p)
    res[1] = du[1]       # ψ'(0)=0
    res[2] = u[2]       # ϕ(0)=0
    # res[3] = u[1] - 1   # ψ(0) = 1: extra constraint to fix unknown parameter Q
end

function bcb!(res, du, u, p)
    res[1] = u[1]       # ψ(L)=0
    res[2] = u[2]       # ϕ(L)=0
end

# Initial guess for solution profile
function initial_guess(p, t)
    [exp(-t^2); t * exp(-t^2)]
end

bvp = TwoPointSecondOrderBVProblem(
    tearing!,
    (bca!, bcb!),
    initial_guess,
    tspan,
    [Q_guess],
    bcresid_prototype=(zeros(3), zeros(2)),
    fit_parameters=true
)

# Solve
sol = solve(bvp, MIRKN4(), dt=0.01)

# print the estimated value of γ which satisfies the BCs
Q = sol.prob.p[1]

plot(sol, idxs=(0, 1), label=L"ψ(t)", continuity=:right)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:topright, continuity=:right)