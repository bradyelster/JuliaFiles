# sho_bvp_n2.jl
# Quantum harmonic oscillator solved as a boundary-value eigenproblem
# using BoundaryValueDiffEq.jl

using BoundaryValueDiffEq
using Plots

# ---- Physical / numerical parameters ----
const L = 8.0 # domain half-width
tspan = (-L, L)

exact2(t) = (1 / (π^(1 / 4) * sqrt(2^(2) * factorial(2)))) * (t^2 - 1) * exp(-0.5 * t^2)

function ode!(du, u, p, t)
    du[1] = u[2]
    du[2] = (t^2 - 2 * p[1]) * u[1]
end

# ---- Boundary conditions ----
function bca!(residual, ua, p, t)
    residual[1] = ua[1] - exact2(-L)
    #residual[2] = u(0.0)[1] - 1
    #residual[3] = u(L)[1] - exact2(L)
end

function bcb!(residual, ub, p, t)
    residual[1] = ub[1] - exact2(L)
end

guess(p, t) = [exp(-t^2 / 2), -2 * t * exp(-t^2 / 2)]  # initial guess is a Gaussian

# initial parameter guess (energy E)
E_guess = 2.50 # 2.5  # analytic for n=2: 2.5

# ---- Define and solve the BVP ----
bvp = BVProblem(ode!, (bca!, bcb!), guess, tspan, [E_guess], bcresid_prototype=(zeros(1), zeros(1)))

@time sol = solve(bvp, MIRK6(), dt=0.05)

println("fitted energy = ", sol.prob.p)

plot(sol)