using BoundaryValueDiffEq, Plots, LaTeXStrings

# --- Physical parameters ---
L = 15.0
tspan = (0.0, L)
f(t) = tanh(t)
ddf(t) = -2 * tanh(t) * sech(t)^2

k = 0.5
#S = 100.0
#ε = S^(-2 / 3)
#Q_guess = 0.26  # initial guess for eigenvalue-like parameter Q

# --- Differential equations ---
function outer!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    du[1] = ψ1
    du[2] =
        du[3] = (k^2 - (2 * sech(t)^2)) * ψ
end

# --- Boundary conditions ---
# Half-domain version assumes ψ even, ϕ odd
function bca!(res, u, p)
    ψ, ψ1 = u
    res[1] = ψ - 1    # ψ'(0)=0 (even)
end

function bcb!(res, u, p)
    ψ, ψ1 = u
    res[1] = ψ - (exp(-k * 15) * (1 + k^(-1) * f(15)))   # ψ(L)=0 (decay)
end

# --- Initial guess ---
function initial_guess(p, t)
    ψ0 = exp(-k * t) * (1 + k^(-1) * f(t))
    ψ1 = (exp(-k * t) * (sech(t)^2 - k * (k + f(t)))) / k
    return [ψ0, ψ1]
end

bvp = TwoPointBVProblem(
    outer!, (bca!, bcb!), initial_guess, tspan;
    bcresid_prototype=(zeros(1), zeros(1))
)

sol = solve(bvp, MIRK4(), dt=0.01, verbose=true, reltol=1e-8, abstol=1e-9)

# --- Plot results ---
plot(sol, idxs=(0, 1), lw=2)