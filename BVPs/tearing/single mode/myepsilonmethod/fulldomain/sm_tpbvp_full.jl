# Tearing Mode - Numerical Solution, 10/01/2025
# Single Helicity Mode, BVProblem

using BoundaryValueDiffEq, Plots, LaTeXStrings

L = 15.0
tspan = (-L, L)
f(t) = tanh(t)
ddf(t) = -2 * tanh(t) * sech(t)^2

k = 0.5
S = 100
ε = S^(-2 / 3)
Q_guess = 0.2607 # 0.18148777711494185 # for S = 1000, # 0.2607 for S=100 # Q -> 1 is the NCF limit where γ = S^(-1/3) Q, CF limit is Q << 1

function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    Q = p[1]
    du[1] = ψ1
    du[2] = ϕ1
    du[3] = (Q / ε + k^2) * ψ - (f(t) / ε) * ϕ
    du[4] = (ddf(t) / (Q * ε) - f(t) / ε^2) * ψ + (f(t)^2 / (Q * ε^2) + k^2) * ϕ
end

function bc!(res, u, p, t)
    res[1] = u(-L)[1]       # ψ(-L)=0
    res[2] = u(-L)[2]       # ϕ(-L)=0
    res[3] = u(L)[1]       # ψ(L)=0
    res[4] = u(L)[2]       # ϕ(L)=0
    res[5] = u(0.0)[1] - 1   # ψ(0) = 1: extra constraint to fix unknown parameter Q
end

# Initial guess for solution profile
function initial_guess(p, t)
    [exp(-t^2), t * exp(-t^2), -2 * t * exp(-t^2), (1 - 2 * t^2) * exp(-t^2)]
end

bvp = BVProblem(
    BVPFunction(tearing!, bc!; bcresid_prototype=zeros(5)),
    initial_guess, tspan, [Q_guess], fit_parameters=true)

@time sol = solve(bvp, MIRK4(), dt=0.01,
    saveat=0.1,
    adaptive=true,
    progress=true,
    verbose=true
)

# print the estimated value of γ which satisfies the BCs
Q = sol.prob.p[1]

plot(sol, idxs=(0, 1), label=L"ψ(t)", continuity=:right)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:topright, continuity=:right)