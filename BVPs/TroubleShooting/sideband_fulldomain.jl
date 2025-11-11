# Tearing Mode Numerical Solution
# Full Domain (-L, L) w/ B&L 2018's Normalizations

using BoundaryValueDiffEq, Plots, LaTeXStrings

L = 15.0
tspan = (-L, L)
# Define the functions to plot
#f_mplus(t) = kp * (1 + 0.1*f(t)) - kn * Bz
#f_mminus(t) = km * (1 + 0.1*f(t)) - kn * Bz
#f_m(t) = k * (1 + 0.1*f(t)) - kn * Bz

@inline f(t) = 0.1 * tanh(t)
@inline ddf(t) = -2 * 0.1 * tanh(t) * sech(t)^2

# mode numbers of interest
m = 2
n = 1

# device parameters
R0 = 1

# poloidal wave numbers 
ky = 0.25
km = (m - 1) * ky
k = (m) * ky
kp = (m + 1) * ky

# toroidal wave number
kz = 1 / R0
kn = n * kz

# combined wave numbers
K2 = k^2 # + kn^2
Km2 = km^2 # + kn^2
Kp2 = kp^2 # + kn^2

# Lundquist number
S = 10

function sideband!(du, u, p, t)
    ψ, ψ1, ϕ, ϕ1, x, γ = u
    du[1] = ψ1
    du[2] = S * γ * ψ - S * k * f(t) * ϕ + K2 * ψ
    du[3] = ϕ1
    du[4] = -((k + ky) * f(t) + ky) * ψ / γ + ((k + ky) * f(t)) * Kp2 * ψ / γ + ((k + ky) * ddf(t) + (ky * ddf(t) / f(t))) * ψ / γ + Kp2 * ϕ
    du[5] = abs2(ψ)
    du[6] = 0
end

# initial state vector at t=-L, informed from half-domain solution
u0 = [1e-5, 1e-5, 0.001, 0.001, 0.0, 0.0661]

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

bvp = TwoPointBVProblem(sideband!, (bca!, bcb!), u0, tspan, bcresid_prototype=(zeros(3), zeros(3)))

@time sol = solve(bvp, MIRK4(), dt=0.05)

γ_found = round(sol.u[1][6], digits=4)

plot(sol, idxs=(0, 1), label=L"ψ(x)", lw=2)
plot!(sol, idxs=(0, 3), label=L"φ(x)", xlabel="x", legend=:best, lw=2, title="γ = $γ_found")
# plot the (reduced) mesh used by solver
scatter!(sol.t[1:4:end], zeros(length(sol.t[1:4:end])), markershape=:vline, color="lightgray", label="mesh")
plot!(sol.t, f.(sol.t), label="mag. field profile")

#=
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
=#