# Tearing Mode Numerical Solution
# Full Domain (-L, L) w/ Fitzpatrick's Normalizations

using BoundaryValueDiffEq, Plots, LaTeXStrings

const L = 20.0
tspan = (-L, L)
@inline f(t) = tanh(t)
@inline ddf(t) = -2 * tanh(t) * sech(t)^2

const k = 0.5
const S = 50
γ_guess = 0.013

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

# initial state vector at t=-L, informed from half-domain solution
u0 = [3.303002113950678e-28,
    -3.154457544752521e-28,
    -0.0014512794069637734,
    -0.0014512794069642813]

function bca!(res, u, p)
    res[1] = u[1] # - ψouter(-L) # ψ(-L) ≈ 0
    res[2] = u[2] # - ϕouter(-L) # ϕ(-L) ≈ 0
    res[3] = u[4]
end

function bcb!(res, u, p)
    res[1] = u[1] # - ψouter(L) # ψ(L) ≈ 0
    res[2] = u[2] # - ϕouter(L) # ϕ(L) ≈ 0
end

bvp = TwoPointBVProblem(tearing!, (bca!, bcb!), u0, tspan, [γ_guess], bcresid_prototype=(zeros(3), zeros(2)), fit_parameters=true)

@time sol = solve(bvp, MIRK6(), dt=0.01,
    reltol=1e-7,
    abstol=1e-7,
    adaptive=true,
    progress=true,
    verbose=true
)

println("γ = ", sol.prob.p[1])

plot(sol, idxs=(0, 1), label=L"ψ(t)", lw=2)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:best, lw=2)

#=
# provide analytic jacobian
function tearing_jac!(J, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    γ = p[1]
    J[1, 1] = 0
    J[1, 2] = 0
    J[1, 3] = 1
    J[1, 4] = 0
    J[2, 1] = 0
    J[2, 2] = 0
    J[2, 3] = 0
    J[2, 4] = 1
    J[3, 1] = S * γ + k^2
    J[3, 2] = -S * γ * f(t)
    J[3, 3] = 0
    J[3, 4] = 0
    J[4, 1] = ((ddf(t) - f(t) * S * γ) / γ^2)
    J[4, 2] = ((k^2 * γ^2 + S * γ * f(t)^2) / γ^2)
    J[4, 3] = 0
    J[4, 4] = 0
    nothing
end

# provide parameter jacobian pJ (returns df/dp for p = [γ])
function paramjac!(pJ, u, p, t)
    γ = p[1]
    ψ, ϕ, ψ1, ϕ1 = u
    # derivatives of A and B w.r.t γ
    dA = -2 * ddf(t) / γ^3 + S * f(t) / γ^2
    dB = -S * f(t)^2 / γ^2
    # pJ is 4 x nparams (here 4x1)
    pJ[1, 1] = 0.0
    pJ[2, 1] = 0.0
    pJ[3, 1] = S * ψ - S * f(t) * ϕ
    pJ[4, 1] = dA * ψ + dB * ϕ
    nothing
end
=#