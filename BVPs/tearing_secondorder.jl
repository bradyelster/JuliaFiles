using BoundaryValueDiffEq, Plots, LaTeXStrings

L = 15.0
tspan = (0.0, L)
F(t) = tanh(t)
Fpp(t) = -2 * tanh(t) * sech(t)^2

k = 0.5
S = 100
γ_guess = S^(-2 / 5)

function tearing!(ddu, u, p, t)
    ψ, φ = u
    γ = p[1]
    ddu[1] = (S * γ + k^2) * ψ - S * γ * F(t) * φ
    ddu[2] = (Fpp(t) / γ^2 - F(t) * S / γ) * ψ + (S * F(t)^2 / γ + k^2) * φ
end

function bc!(res, u, p)
    res[1] = du[1][1]       # ψ′(0) = 0 
    res[2] = u[1][2]       # φ(0) = 0
    res[3] = u[end][1]     # ψ(L) = 0
    res[4] = u[end][2]     # φ(L) = 0
end

function initial_guess(p, t)
    [exp(-t^2), t * exp(-t^2)]
end

u0 = [1.0, 0.0]

bvp = SecondOrderBVProblem(
    tearing!,
    bc!,
    initial_guess,
    tspan,
    [γ_guess]
)

sol = solve(bvp, MIRKN4(), dt=0.05)

γ = sol.prob.p[1]
println("γ: ", γ)

plot(sol, idxs=(0, 1), label=L"\psi(t)", continuity=:right)
plot!(sol, idxs=(0, 2), label=L"\phi(t)", xlabel="t", legend=:topright, continuity=:right)