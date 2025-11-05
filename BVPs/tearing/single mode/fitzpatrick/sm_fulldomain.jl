# Tearing Mode Numerical Solution
# Full Domain (-L, L) w/ Fitzpatrick's Normalizations

using BoundaryValueDiffEq, Plots, LaTeXStrings

const L = 15.0
tspan = (-L, L)
@inline f(t) = tanh(t)
@inline ddf(t) = -2 * tanh(t) * sech(t)^2

const k = 0.5
const S = 10

function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1, γ, x = u
    du[1] = ψ1
    du[2] = ϕ1
    du[3] = (S * γ + k^2) * ψ - (S * γ * f(t)) * ϕ
    du[4] = ((ddf(t) - f(t) * S * γ) / γ^2) * ψ + ((k^2 * γ^2 + S * γ * f(t)^2) / γ^2) * ϕ
    du[5] = 0
    du[6] = abs2(ψ) # x(t) = |ψ|^2
end

# initial state vector at t=-L, informed from half-domain solution
u0 = [1e-5, 1e-5, 0.001, 0.001, 0.1, 0.0]

function bca!(res, u, p)
    res[1] = u[1] # - ψouter(-L) # ψ(-L) ≈ 0
    res[2] = u[2] # - ϕouter(-L) # ϕ(-L) ≈ 0
    res[3] = u[6]
end

function bcb!(res, u, p)
    res[1] = u[1] # - ψouter(-L) # ψ(-L) ≈ 0
    res[2] = u[2] #- ϕouter(-L) # ϕ(-L) ≈ 0
    res[3] = u[6] - 1
end

bvp = TwoPointBVProblem(tearing!, (bca!, bcb!), u0, tspan, bcresid_prototype=(zeros(3), zeros(3)))

@time sol = solve(bvp, MIRK4(), dt=0.1, tstops=[0.0])

γ_found = round(sol.u[1][5], digits=4)

plot(sol, idxs=(0, 1), label=L"ψ(t)", lw=2)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:best, lw=2, title="γ = $γ_found")
# plot the (reduced) mesh used by solver
scatter!(sol.t[1:4:end], zeros(length(sol.t[1:4:end])), markershape=:vline, color="lightgray", label="mesh")

# plot more physically relevant quantities 
x = sol.t                # radial coordinate (your independent var)
ψ = [u[1] for u in sol.u]  # magnetic flux eigenfunction
ϕ = [u[2] for u in sol.u]  # velocity streamfunction
γ = sol.u[1][5]             # eigenvalue

# set up a symmetric domain around x=0, y=0
Nx, Ny = 400, 400
xg = range(-maximum(abs.(x)), stop=maximum(abs.(x)), length=Nx)
yg = range(-3π, 3π, length=Ny)

# interpolate ψ(x) and ϕ(x) onto uniform grid
using Interpolations
ψ_itp = interpolate((x,), ψ, Gridded(Linear()))
ϕ_itp = interpolate((x,), ϕ, Gridded(Linear()))

# create 2D real physical fields
Ψ = [real(ψ_itp(xi) * exp(im * k * yj)) for xi in xg, yj in yg]
Φ = [real(ϕ_itp(xi) * exp(im * k * yj)) for xi in xg, yj in yg]

# magnetic streamlines 
contourf(xg, yg, Ψ', aspect_ratio=1, c=:viridis, title="Magnetic Flux Streamlines")

plot!(sol,
    idxs=(0, 1),
    label=L"ψ(x)",
    lw=2,
    ylims=(yg[begin], yg[end]),
    xlims=(-15, 15),
    xlabel="x",
    ylabel="ψ(x)"
)

# velocity streamlines 
contourf(xg, yg, Φ', aspect_ratio=1, c=:plasma, title="Velocity Streamlines")
plot!(sol, idxs=(0, 2), label=L"φ(x)",
    legend=:best,
    lw=2,
    ylims=(yg[begin], yg[end]),
    #xlims=(-5, 5),
    #color="black",
    xlabel="x",
    ylabel="y"
)
