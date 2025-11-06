using BoundaryValueDiffEq, Plots, LaTeXStrings
# using LinearAlgebra # for mass matrix

L = 15.0
tspan = (-L, L)
f(t) = tanh(t)
g(t) = -2 * sech(t)^2

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
Km2 = km^2 + kn^2
K2 = k^2 + kn^2
Kp2 = kp^2 + kn^2

# Zeta-components
ζz = 1
ζy = 1

# Lundquist number
S = 10

# Alfvénic Mach number
# M = 1 / 2 #S^(-1 / 2)

# model the m+1 sideband mode 
function sideband!(du, u, p, t)
    ϕp, ψp, ψp1, ϕp1, x, γ = u
    du[1] = ϕp1
    du[2] = ψp1
    du[3] = S * ((Kp2 / S) * ψp + γ * ψp - ky * (1 + (1 + m) * f(t)) * ϕp)
    du[4] = (1 / γ) * (-Kp2 * γ * ϕp - (Kp2 * ψp + g(t) * ψp - S * ((Kp2 / S) * ψp + γ * ψp - ky * ((1 + 1 + m) * f(t)) * ϕp)) * ky * (1 + (1 + m) * f(t)))
    du[5] = abs2(ψp) # x(t) = |ψ|^2
    du[6] = 0
end

# initial state vector at t=-L, informed from half-domain solution
u0 = [1e-5, 1e-5, 0.001, 0.001, 0.0, 0.123]

# ϕp, ψp, ψp1, ϕp1, ϕp2, x, γ = u
function bca!(res, u, p)
    res[1] = u[1]
    res[2] = u[2]
    res[3] = u[5]
end

function bcb!(res, u, p)
    res[1] = u[1]
    res[2] = u[2]
    res[3] = u[5] - 1
end

fun = BVPFunction(sideband!, (bca!, bcb!); bcresid_prototype=(zeros(3), zeros(3)), twopoint=Val(true))
prob = BVProblem(fun, u0, tspan)

@time sol = solve(prob, MIRK6(), dt=0.08, tstops=[0.0])

plot(sol, idxs=(0, 1), label=L"\phi_{m+1}")
plot!(sol, idxs=(0, 2), label=L"\psi_{m+1}")

tvals = LinRange(-L, L, 1000)
plot!(tvals, kp .* (1 .+ f.(tvals)) .- kn * 1, label=L"\vec{k}_{m+1} \cdot \vec{B}_{0}")
# plot!(tvals, km .* (1 .+ f.(tvals)) .- kn * 1, label="m-1")
# plot!(tvals, k .* f.(tvals), label="m")