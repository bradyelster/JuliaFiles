using BoundaryValueDiffEq, Plots, LaTeXStrings

# Domain
L = 20.0
tspan = (0.0, L)

# Equilibrium profile
F(t) = tanh(t)
Fpp(t) = -2 * tanh(t) * sech(t)^2

# Parameters
k = 0.5
S = 200
γ_guess = S^(-1/3)

# System of equations
function tearing!(du, u, p, t)
    ψ, φ, ψp, φp = u
    γ = p[1]

    # First-order evolution
    du[1] = ψp
    du[2] = φp

    # ψ'' equation
    ψpp = (S * γ) * (ψ - F(t) * φ) + k^2 * ψ
    du[3] = ψpp

    # φ'' equation
    du[4] = (-F(t)/γ^2) * (ψpp - k^2*ψ - (Fpp(t)/F(t))*ψ) + k^2*φ
end

# Boundary conditions
function bca!(res, u, p)
    res[1] = u[3]       # ψ'(0) = 0
    res[2] = u[2]       # φ(0)  = 0
end

function bcb!(res, u, p)
    res[1] = u[1]       # ψ(L) = 0
    res[2] = u[2]       # φ(L) = 0
end

# Initial guess for solution profile
function initial_guess(p, t)
    [exp(-t^2); t * exp(-t^2); -2 * t * exp(-t^2); (1 - 2 * t^2) * exp(-t^2)]
end

# Define BVP
bvp = TwoPointBVProblem(
    tearing!,
    (bca!, bcb!),
    initial_guess,
    tspan,
    [γ_guess],
    bcresid_prototype=(zeros(2), zeros(2)),
    fit_parameters=true
)

# Solve
sol = solve(bvp, MIRK4(), dt=0.05)

# Print estimated γ
println("γ: ", sol.prob.p[1])

# Plot solutions
plot(sol, idxs=(0, 1), label=L"\psi(t)", continuity=:right)
plot!(sol, idxs=(0, 2), label=L"\phi(t)", xlabel="t", legend=:topright, continuity=:right)
