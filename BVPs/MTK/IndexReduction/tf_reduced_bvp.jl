using BoundaryValueDiffEq, LinearAlgebra, Plots

const L = 20.0
tspan = (-L, L) # full-domain

f(t) = tanh(t)
df(t) = 1 - tanh(t)^2 # = sech(t)^2
ddfoverf(t) = -2 * sech(t)^2
ddfdfoverf2(t) = -2 * csch(t) * sech(t)^3
dddfoverf(t) = 4 * tanh(t) * sech(t)^2 - 2 * csch(t) * sech(t)^3

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
S = 100

# Alfvénic Mach number
M = S^(-1 / 2)

function tftearing!(du, u, p, t)
    # parameter to fit: γ
    γ = p[1]

    # state vector elements
    ψp, φp, ψ,
    ψp1, ψm, ψ1,
    ψm1, φm, φ,
    φ1, φ2, φp1,
    φm1, φ3, φm3 = u

    # helper functions/"observables"
    ψp3 = ((Kp2 / S) * ψp1 + γ * ψp1 - ky * (1 + m) * df(t) * φp - ky * (1 + (1 + m) * f(t)) * φp1 + 0.5 * M * S * ((K2 / S) * ψ + 0.5 * M * (ψp1 + ψm1) + γ * ψ - k * f(t) * φ)) * S
    ψm3 = ((Km2 / S) * ψm1 + γ * ψm1 + ky * (1 - m) * df(t) * φm - ky * (-1 + (-1 + m) * f(t)) * φm1 + 0.5 * M * S * ((K2 / S) * ψ + 0.5 * M * (ψp1 + ψm1) + γ * ψ - k * f(t) * φ)) * S
    φp3 = f(t) * (1 / 0.5 * M) * ((-φ2 + K2 * φ) * γ - 0.5 * M * (φm3 - (Kp2 - ky^2) * φm1 - (Kp2 - ky^2) * φp1) - (-ddfoverf(t) * ψ - K2 * ψ + S * ((K2 / S) * ψ + 0.5 * M * (ψp1 + ψm1) + γ * ψ - k * f(t) * φ)) * k * f(t))

    # ̇v = A v setup
    du[1] = ψp1
    du[2] = φp1
    du[3] = ψ1
    du[4] = Kp2 * ψp + 0.5 * M * S * ψ1 + S * γ * ψp - S * ky * (1 + (1 + m) * f(t)) * φp
    du[5] = ψm1
    du[6] = K2 * ψ + 0.5 * M * S * (ψp1 + ψm1) + S * γ * ψ - S * k * f(t) * φ
    du[7] = Km2 * ψm + 0.5 * M * S * ψ1 + S * γ * ψm - ky * (-1 + (-1 + m) * f(t)) * φm
    du[8] = φm1
    du[9] = φ1
    du[10] = φ2
    du[11] = φ3
    du[12] = f(t) * (-1 / γ) * (-Kp2 * γ * φ1 + 0.5 * M * (φ3 + (-Kp2 + ζz^2 * ky^2) * φ1) + (-ddfoverf(t) * ψp - Kp2 * ψp + S * (Kp2 * ψp / S + 0.5 * M * ψ1 + ψp * γ - ky * (1 + (1 + m) * f(t)) * φp)) * ky * (1 + (1 + m) * f(t)))
    du[13] = f(t) * (-1 / γ) * (-Km2 * γ * φ1 + 0.5 * M * (φ3 + (-Km2 + ζz^2 * ky^2) * φ1) + (-ddfoverf(t) * ψm - Km2 * ψm + S * (Km2 * ψm / S + 0.5 * M * ψ1 + ψm * γ - ky * (-1 + (-1 + m) * f(t)) * φm)) * ky * (-1 + (-1 + m) * f(t)))
    du[14] = f(t)^2 * (-1 / 0.5 * M) * (-(φm3 + Km2 * φm1) * γ - (dddfoverf(t) * ψm + ddfoverf(t) * ψm1 - ddfdfoverf2(t) * ψm - ψm3 + Km2 * ψm1) * ky * (-1 + (-1 + m) * f(t)) + 0.5 * M * (-Km2 + ζz^2 * ky^2) * φ2 + (ddfoverf(t) * ψm + Km2 * ψm + S * (-Km2 * ψm / S - 0.5 * M * ψ1 - γ * ψm + ky * (-1 + (-1 + m) * f(t)) * φm)) * ky * (1 - m) * df(t))
    du[15] = f(t)^2 * ((-φp3 + Kp2 * φp1) * γ + ((-ddfdfoverf2(t) * ψp) - ψp3 + (dddfoverf(t) * ψp + ddfoverf(t) * ψp1) + Kp2 * ψp1) * ky * (1 + (1 + m) * f(t)) - 0.5 * M * (-(-φm3 + Km2 * φm1) * γ - ((dddfoverf(t) * ψm + ddfoverf(t) * ψm1) + (-ψm * ddfdfoverf2(t) - ψm3 + Km2 * ψm1) * ky * (-1 + (-1 + m) * f(t)) + 0.5 * M * (-Km2 + (ky^2) * (ζz^2)) * φ2 + (ψm * ddfoverf(t) + Km2 * ψm + S * ((-Km2 * ψm) / S - 0.5 * M * ψ1 - ψm * γ + ky * (-1 + (-1 + m) * f(t)) * φm)) * ky * (1 - m) * df(t)) / (0.5 * M) + (-Kp2 + (ky^2) * (ζz^2)) * φ2) + (ddfoverf(t) * ψp + Kp2 * ψp + S * ((-Kp2 * ψp) / S - 0.5 * M * ψ1 - ψp * γ + ky * (1 + (1 + m) * f(t)) * φp)) * ky * (1 + m) * df(t))
end

mat = Matrix{Float64}(I, 15, 15)  # creates a 15×15 identity matrix
mat[15, 15] = 0                   # set the last element to 0

function bc!(res, u, p, t)
    # ψp, φp, ψ,
    # ψp1, ψm, ψ1,
    # ψm1, φm, φ,
    # φ1, φ2, φp1,
    # φm1, φ3, φm3 = u

    res[1] = u(-L)[3] - 1e-5     # ψm(L)=0,    Dirichlet on left boundary
    res[2] = u(-L)[5] - 1e-5    # ψm-1(L)=0,  Dirichlet on left boundary
    res[3] = u(-L)[1] - 1e-5   # ψm+1(L)=0,  Dirichlet on left boundary
    res[4] = u(-L)[9] - 1e-5    # ϕm(L)=0,    Dirichlet on left boundary
    res[5] = u(-L)[8] - 1e-5    # ϕm-1(L)=0,  Dirichlet on left boundary
    res[6] = u(-L)[2] - 1e-5    # ϕm+1(L)=0,  Dirichlet on left boundary
    res[7] = u(L)[3] - 1e-5    # ψm(L)=0,    Dirichlet on left boundary
    res[8] = u(L)[5] - 1e-5    # ψm-1(L)=0,  Dirichlet on left boundary
    res[9] = u(L)[1] - 1e-5   # ψm+1(L)=0,  Dirichlet on left boundary
    res[10] = u(L)[9] - 1e-5   # ϕm(L)=0,    Dirichlet on left boundary
    res[11] = u(L)[8] - 1e-5   # ϕm-1(L)=0,  Dirichlet on left boundary
    res[12] = u(L)[2] - 1e-5   # ϕm+1(L)=0,  Dirichlet on left boundary
    res[13] = u(0.0)[6] - 1e-5   # ψm'(0) = 0,  even
    res[14] = u(0.0)[9] - 1e-5   # ϕm(L)=0
    #res[15] = u(0.0)[3] - 1  # normalization on ψm 
end

u0 = [0.0, 0.0, 0.0,
    0.01, 0.0, 0.01,
    0.0, 0.0, 0.0,
    0.01, 0.0, 0.01,
    0.01, 0.0]

γ_guess = 0.05
fun = BVPFunction(tftearing!, bc!; mass_matrix=mat, bcresid_prototype=zeros(14))
prob = BVProblem(fun, u0, tspan, [γ_guess], fit_parameters=false)

sol = solve(prob, MIRK4(), dt=0.1)

plot(sol, idxs=(0, 3), label="ψ")
plot!(sol, idxs=(0, 9), label="φ")