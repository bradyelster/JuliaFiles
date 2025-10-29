# Tearing Mode Numerical Solution
# Full Domain (-L, L) w/ Fitzpatrick's Normalizations

using BoundaryValueDiffEq, Plots, LaTeXStrings

const L = 15.0
tspan = (-L, L)
@inline f(t) = tanh(t)
@inline ddf(t) = -2 * tanh(t) * sech(t)^2

const k = 0.5
const S = 10
γ_guess = 0.1

function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    γ = p[1]
    du[1] = ψ1
    du[2] = ϕ1
    du[3] = (S * γ + k^2) * ψ - (S * γ * f(t)) * ϕ
    du[4] = ((ddf(t) - f(t) * S * γ) / γ^2) * ψ + ((k^2 * γ^2 + S * γ * f(t)^2) / γ^2) * ϕ
end

# initial state vector at t=-L, informed from half-domain solution
u0 = [1e-5, 1e-5, 0.001, 0.001]

function bc!(res, u, p, t)
    res[1] = u(-L)[1] # - ψouter(-L) # ψ(-L) ≈ 0
    res[2] = u(-L)[2] # - ϕouter(-L) # ϕ(-L) ≈ 0
    #res[3] = u(0.0)[1] - 1
    res[4] = u(L)[1] # - ψouter(L) # ψ(L) ≈ 0
    res[5] = u(L)[2] # - ϕouter(L) # ϕ(L) ≈ 0
end

bvp = BVProblem(tearing!, bc!, u0, tspan, [γ_guess], fit_parameters=true)

@time sol = solve(bvp, Ascher4(zeta=[0.0, 0.0, 0.5, 1.0, 1.0]), dt=0.05)

println("γ = ", sol.prob.p[1])

plot(sol, idxs=(0, 1), label=L"ψ(t)", lw=2)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:best, lw=2)