# Tearing Mode Numerical Solution
# Full Domain (-L, L) w/ Fitzpatrick's Normalizations

using BoundaryValueDiffEq, Plots, LaTeXStrings

const L = 15.0
tspan = (0.0, L)
@inline f(t) = tanh(t)
@inline ddf(t) = -2 * tanh(t) * sech(t)^2

const k = 0.5
const S = 100
γ_guess = 0.05753900476044987

# provide outer solutions for easier boundary value handling
ψouter(t) = exp(-k * abs(t)) * (1 + f(abs(t)) / k)
ϕouter(t) = ψouter(t) / f(t)
ψouter_prime(t) = -k * sign(t) * exp(-k * abs(t)) * (1 + f(abs(t)) / k) + exp(-k * abs(t)) * (ddf(abs(t)) / k)
ϕouter_prime(t) = ψouter_prime(t) / f(t) - ψouter(t) * ddf(t) / (f(t)^2)

function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    γ = p[1]
    du[1] = ψ1
    du[2] = ϕ1
    du[3] = (S * γ + k^2) * ψ - (S * γ * f(t)) * ϕ
    du[4] = ((ddf(t) - f(t) * S * γ) / γ^2) * ψ + ((k^2 * γ^2 + S * γ * f(t)^2) / γ^2) * ϕ
end

# initial state vector at t=0
u0 = [1.0, 0.0, 0.01, 1.0]

function bca!(res, u, p)
    res[1] = u[3] # ψ'(0) = 0
    res[2] = u[2] # ϕ(0) = 0
    res[3] = u[1] - 1 # normalization : ψ(0) = 1
end

function bcb!(res, u, p)
    res[1] = u[1] # ψ(L) 0 0
    res[2] = u[2] # ϕ(L) 0 0
end

bvp = TwoPointBVProblem(tearing!, (bca!, bcb!), u0, tspan, [γ_guess], bcresid_prototype=(zeros(3), zeros(2)), fit_parameters=true)

@time sol = solve(bvp, MIRK4(), dt=0.01,
    adaptive=true,
    progress=true,
    verbose=true
)

println("γ = ", sol.prob.p[1])

plot(sol, idxs=(0, 1), label=L"ψ(t)", lw=2)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:best, lw=2)

sol.u[1]

# Extract time and solution arrays
t_half = sol.t
u_half = hcat(sol.u...)  # columns correspond to time points
ψ_half = u_half[1, :]
ϕ_half = u_half[2, :]
ψp_half = u_half[3, :]  # ψ'(t)
ϕp_half = u_half[4, :]  # ϕ'(t)


# --- Extend to full domain ---
# (Don't double count t=0)
t_full = [-reverse(t_half[2:end]); t_half]
ψ_full = [reverse(ψ_half[2:end]); ψ_half]        # even extension
ϕ_full = [-reverse(ϕ_half[2:end]); ϕ_half]       # odd extension
ψp_full = [-reverse(ψp_half[2:end]); ψp_half]       # ψ even → ψ' odd
ϕp_full = [reverse(ϕp_half[2:end]); ϕp_half]        # ϕ odd  → ϕ' even

# --- Extract final state vector at t = L ---
final_u = sol.u[end]

# --- plot full-domain extensions ---
# plot(t_full, ψ_full, label=L"\psi(t)", lw=2, title="Fitzpatrick: Tearing Mode Solution")
# plot!(t_full, ϕ_full, label=L"\phi(t)", lw=2, xlabel="x", legend=:best)

# --- START FULL DOMAIN SOLUTION ATTEMPT --- #
tspan2 = (-L, L)
γ_guess2 = sol.prob.p[1] # use guess from half-domain solution

# guess = [ψ_full, ϕ_full, ψp_full, ϕp_full]

function bca2!(res, u, p)
    res[1] = u[1]
    res[2] = u[2]
end

function bcb2!(res, u, p)
    res[1] = u[1] # - ψouter(L) # ψ(L) ≈ 0
    res[2] = u[2] # - ϕouter(L) # ϕ(L) ≈ 0
end

bvp2 = TwoPointBVProblem(
    BVPFunction(tearing!, (bca2!, bcb2!); twopoint=Val(true), bcresid_prototype=(zeros(2), zeros(2))),
    final_u, tspan2, [γ_guess2], fit_parameters=true
)

@time sol2 = solve(bvp2, MIRK6(), dt=0.01,
    reltol=1e-7,
    abstol=1e-7,
    adaptive=true,
    progress=true,
    verbose=true
)

println("γ = ", sol2.prob.p[1])
sol2.prob.p[1] - sol.prob.p[1]

plot(sol2, idxs=(0, 1), label=L"ψ(t)", lw=2)
plot!(sol2, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:best, lw=2)


# --- END FULL DOMAIN SOLUTION ATTEMPT --- #