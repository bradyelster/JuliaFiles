using BoundaryValueDiffEq, Plots, LaTeXStrings

L = 15.0
tspan = (-L, L)
f(t) = tanh(t)

k = 0.5
S = 10

function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1, x, γ = u

    du[1] = ψ1
    du[2] = ϕ1
    du[3] = (-S * k * f(t)) * ϕ + (S * γ - k^2) * ψ
    du[4] = ((k^2 * S * f(t)^2) / γ) * ϕ - (k * f(t) / γ) * (-2 * k^2 + S * γ + 2 * (sech(t)^2)) * ψ
    du[5] = abs2(ψ)
    du[6] = 0
end

# initial state vector at t=-L, informed from half-domain solution
u0 = [1e-5, 1e-5, 0.001, 0.001, 0.0, 0.123]

function bca!(res, u, p)
    res[1] = u[1] # - ψouter(-L) # ψ(-L) ≈ 0
    res[2] = u[2] # - ϕouter(-L) # ϕ(-L) ≈ 0
    res[3] = u[5]
end

function bcb!(res, u, p)
    res[1] = u[1] # - ψouter(-L) # ψ(-L) ≈ 0
    res[2] = u[2] #- ϕouter(-L) # ϕ(-L) ≈ 0
    res[3] = u[5] - 1
end

bvp = TwoPointBVProblem(tearing!, (bca!, bcb!), u0, tspan, bcresid_prototype=(zeros(3), zeros(3)))

@time sol = solve(bvp, MIRK4(), dt=0.1, tstops=[0.0])

γ_found = round(sol.u[1][6], digits=4)

plot(sol, idxs=(0, 1), label=L"ψ(t)", lw=2)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:best, lw=2, title="γ = $γ_found")
scatter!(sol.t[1:4:end], zeros(length(sol.t[1:4:end])), markershape=:vline, color="lightgray", label="mesh")

#savefig("bl2018.png")

# flip signs of ψ and φ
#plot(sol, x -> -x)