# sho_bvp_n2.jl
# Quantum harmonic oscillator solved as a boundary-value eigenproblem
# using BoundaryValueDiffEq.jl

using BoundaryValueDiffEq
using Plots

# ---- Physical / numerical parameters ----
const L = 8.0 # domain half-width
tspan = (-L, L)

function ode!(du, u, p, t)
    du[1] = u[2]
    du[2] = (t^2 - 2 * p[1]) * u[1]
end

# ---- Boundary conditions ----
function bc!(res, u, p, t)
    res[1] = u(-L)[1] # ψ(-L) ≈ 0 decay
    res[2] = u(0.0)[1] - 1 # normalization condition
    res[3] = u(L)[1] # ψ(L) ≈ 0 decay
end

# ---- Initial guess for the solution ----
guess(p, t) = [exp(-t^2 / 2); -t * exp(-t^2 / 2)]  # initial guess is a Gaussian

# initial parameter guess (energy E)
E_guess = 2.6 # 2.5  # analytic for n=2: 2.5

# ---- Define and solve the BVP ----
bvp = BVProblem(
    BVPFunction(ode!, bc!; bcresid_prototype=zeros(3)),
    guess, tspan, [E_guess], fit_parameters=true)

@time sol = solve(bvp, MIRK6(), dt=0.5, adaptive=true)

println("fitted energy = ", sol.prob.p)

plot(sol)