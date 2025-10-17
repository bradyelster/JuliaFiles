using BoundaryValueDiffEq, Plots

const L = 1
tspan = (0, L)
S = 100
k = 0.5
Δprime(ξ) = 2*(1-ξ^2)/ξ
γ_guess = 0.057

function layer!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    γ = p[1]
    du[1] = ψ1
    du[2] = ϕ1
    du[3] = S*γ*ψ - t*S*γ*ϕ
    du[4] = -t*S*ψ / γ + t^2*S*ϕ / γ
end

# initial state vector at t=-L
u0 = [1.0, 0.0, 0.01, 1.0]

function bca!(res, u, p)
    res[1] = u[3]
    res[2] = u[2]
    res[3] = u[1] - 1 # normalization : ψ(0) = 1
end

function bcb!(res, u, p)
    res[1] = u[1] - (1 + 0.5*Δprime(k)*abs(L))
    res[2] = u[2] - (1 + 0.5*Δprime(k)*abs(L)) / L
end

bvp = TwoPointBVProblem(layer!, (bca!, bcb!), u0, tspan, [γ_guess], bcresid_prototype=(zeros(3), zeros(2)), fit_parameters=true)

@time sol = solve(bvp, MIRK4(), dt=0.01,
    adaptive=true,
    progress=true,
    verbose=true
)

println("γ = ", sol.prob.p[1])

plot(sol, idxs=(0, 1), label="ψ(t)", lw=2)
plot!(sol, idxs=(0, 2), label="φ(t)", xlabel="t", legend=:best, lw=2)