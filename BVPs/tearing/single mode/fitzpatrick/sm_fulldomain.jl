# Tearing Mode Numerical Solution
# Full Domain (-L, L) w/ Fitzpatrick's Normalizations

using BoundaryValueDiffEq, Plots

L = 15.0
tspan = (-L, L)
@inline f(t) = tanh(t)
@inline ddf(t) = -2 * tanh(t) * sech(t)^2

k = 0.5 # Harris sheet is unstable only for k ∈ (0, 1), stable for k > 1. 
S = 10

function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1, x, γ = u
    du[1] = ψ1
    du[2] = ϕ1
    du[3] = (S * γ + k^2) * ψ - (S * γ * f(t)) * ϕ
    du[4] = ((ddf(t) - f(t) * S * γ) / γ^2) * ψ + ((k^2 * γ^2 + S * γ * f(t)^2) / γ^2) * ϕ
    #du[4] = -f(t)*du[3]/γ^2 + (k^2/γ^2)*ψ + (ddf(t)/γ^2)*ψ
    du[5] = abs2(ψ) # x(t) = |ψ|^2
    du[6] = 0
end

# initial state vector at t=-L, informed from half-domain solution
sol = [1e-5, 1e-5, 0.001, 0.001, 0.0, 0.123]

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

bvp = TwoPointBVProblem(tearing!, (bca!, bcb!), sol, tspan, bcresid_prototype=(zeros(3), zeros(3)))
sol = solve(bvp, MIRK4(), dt=0.1, tstops=[0.0])

γ_found = round(sol.u[1][6], digits=4)

plot(sol, idxs=(0, 1), label="", lw=2, lc="black")
plot!(sol, idxs=(0, 2), label="", xlabel="", lw=2, lc=RGB(34 / 255, 88 / 255, 52 / 255))
scatter!(
    sol.t[1:4:end],
    zeros(length(sol.t[1:4:end])),
    markershape=:vline,
    color="black",
    label="",
    background_color=:transparent,
    dpi=300,
    legend=:false
)

savefig("fitzpatrick_S10_v2.png")

#=
# plot more physically relevant quantities 
x = sol.t                # radial coordinate (your independent var)
ψ = [u[1] for u in sol.u]  # magnetic flux eigenfunction
ϕ = [u[2] for u in sol.u]  # velocity streamfunction
γ = sol.u[1][6]             # eigenvalue

# set up a symmetric domain around x=0, y=0
Nx, Ny = 400, 400
xg = range(-maximum(abs.(x)), stop=maximum(abs.(x)), length=Nx)
yg = range(-3π, 3π, length=Ny)

# interpolate ψ(x) and ϕ(x) onto uniform grid
using Interpolations
ψ_itp = interpolate((x,), ψ, Gridded(Linear()))
ϕ_itp = interpolate((x,), ϕ, Gridded(Linear()))

# create 2D real physical fields
Ψ = [-log(cosh(xi)) + real(ψ_itp(xi) * exp(im * k * yj)) for xi in xg, yj in yg]
Φ = [real(ϕ_itp(xi) * exp(im * k * yj)) for xi in xg, yj in yg]
# 0.1 * (1 / k) * sin(k * yj)
# + real(ϕ_itp(xi) * exp(im * k * yj))
# magnetic streamlines 
contourf(xg, yg, Ψ', aspect_ratio=1, c=:viridis, title="Magnetic Flux Streamlines")

# velocity streamlines 
contourf(xg, yg, Φ', c=:plasma, title="Velocity Streamlines", xlims=(-5, 5))
savefig("pseudocontour.png")
plot!(sol, idxs=(0, 2), label=L"φ(x)",
    legend=:best,
    lw=2,
    ylims=(yg[begin], yg[end]),
    xlims=(-5, 5),
    #color="black",
    xlabel="x",
    ylabel="y"
)

plot!(sol,
    idxs=(0, 1),
    label=L"ψ(x)",
    lw=2,
    ylims=(yg[begin], yg[end]),
    xlims=(-15, 15),
    xlabel="x",
    ylabel="ψ(x)"
)
=#
