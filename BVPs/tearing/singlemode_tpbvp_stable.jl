using BoundaryValueDiffEq, Plots, LaTeXStrings

# --- Physical parameters ---
L = 15.0
tspan = (0.0, L)          # Use half-domain if symmetry known
f(t) = tanh(t) + 1
ddf(t) = -2 * tanh(t) * sech(t)^2

k = 0.5
S = 100.0
ε = S^(-2 / 3)
Q_guess = 0.26  # initial guess for eigenvalue-like parameter Q

# --- Differential equations ---
function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    Q = p[1]
    du[1] = ψ1
    du[2] = ϕ1
    du[3] = (Q / ε + k^2) * ψ - (f(t) / ε) * ϕ
    du[4] = (ddf(t) / (Q * ε) - f(t) / ε^2) * ψ + (f(t)^2 / (Q * ε^2) + k^2) * ϕ
end

# --- Boundary conditions ---
# Half-domain version assumes ψ even, ϕ odd
function bca!(res, u, p)
    ψ, ϕ, ψ1, ϕ1 = u
    # symmetry BCs at x=0
    res[1] = ψ1     # ψ'(0)=0 (even)
    res[2] = ϕ      # ϕ(0)=0 (odd)
end

function bcb!(res, u, p)
    ψ, ϕ, ψ1, ϕ1 = u
    res[1] = ψ      # ψ(L)=0 (decay)
    res[2] = ϕ      # ϕ(L)=0 (decay)
end

# --- Initial guess ---
function initial_guess(p, t)
    ψ0 = exp(-t^2 / 10)
    ϕ0 = t * exp(-t^2 / 10)
    ψ1 = -2t / 10 * exp(-t^2 / 10)
    ϕ1 = exp(-t^2 / 10) - 2t^2 / 10 * exp(-t^2 / 10)
    return [ψ0, ϕ0, ψ1, ϕ1]
end

# --- Normalization condition ---
# We add an extra residual enforcing ψ(0) = 1
function normalization_condition!(res, u, p)
    res[1] = u[1] - 1.0
end

function bca!(res, u, p)
    ψ, ϕ, ψ1, ϕ1 = u
    res[1] = ψ1         # ψ'(0)=0
    res[2] = ϕ          # ϕ(0)=0
    res[3] = ψ - 1.0    # normalization ψ(0)=1
end

function bcb!(res, u, p)
    ψ, ϕ, ψ1, ϕ1 = u
    res[1] = ψ
    res[2] = ϕ
end

bvp = TwoPointBVProblem(
    tearing!, (bca!, bcb!), initial_guess, tspan, [Q_guess];
    bcresid_prototype=(zeros(3), zeros(2)),  # match sizes of bca, bcb
    fit_parameters=true
)

@time sol = solve(bvp, MIRK4(), dt=0.01, verbose=true)

Q = sol.prob.p[1]
println("Computed eigenvalue Q ≈ ", Q)

# --- Plot results ---
plot(sol, idxs=(0, 1), label=L"\hat{\psi}(t)", lw=2)
plot!(sol, idxs=(0, 2), label=L"\hat{\phi}(t)", lw=2, xlabel="t", legend=:topright)