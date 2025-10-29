using LinearAlgebra
using BoundaryValueDiffEq
using Plots

const L = 15.0
tspan = (0.0, L)
@inline f(t) = tanh(t)
@inline ddfoverf(t) = -2 * sech(t)^2
@inline safe_ddfoverf(t) = t < 1e-6 ? -2 + 4 / 3 * t^2 : ddfoverf(t)

# mode numbers of interest
m, n = 2, 1

# device parameters
R0 = 1

# poloidal wave numbers 
ky = 0.25
km = (m - 1) * ky
k = (m) * ky
kp = (m + 1) * ky

# toroidal wave number
kz = 1 / R0
kn = (n) * kz

# combined wave numbers
Km2 = km^2 + kn^2
K2 = k^2 + kn^2
Kp2 = kp^2 + kn^2

# Zeta-components (CHANGE THESE LATER)
ζz = 1
ζy = 1

# Lundquist number
S = 10

# Growth rate
γ_guess = 0.0575

# Alfvénic Mach number
M = S^(-1 / 2)

function tftearing!(du, u, p, t)
    phm0, phm1,
    ph0, ph1,
    php0, php1,
    psm0, psm1,
    ps0, ps1,
    psp0, psp1,
    phm2, ph2, php2 = u
    γ = p[1]
    du[1] = phm1
    du[2] = phm2
    du[3] = ph1
    du[4] = ph2
    du[5] = php1
    du[6] = php2
    du[7] = psm1
    du[8] = (S * γ + Km2) * psm0 + S * 0.5 * M * ps1 - (ky * (f(t) * (m - 1) - 1)) * phm0 # ψm-1''
    du[9] = ps1
    du[10] = (S * γ + K2) * ps0 + S * 0.5 * M * (psp1 + psm1) - (S * k * f(t)) * ph0 # ψm''
    du[11] = psp1
    du[12] = (S * γ + Kp2) * psp0 + S * 0.5 * M * ps1 - (ky * (f(t) * (m + 1) + 1)) * php0 # ψm+1''
    du[13] = (-0.5 * ky / M) * (f(t) * (m - 1) - 1) * (du[8] - (Km2 + safe_ddfoverf(t)) * psm0) - (0.5 * γ / M) * (phm2 - Km2 * phm0) + (K2 - ky^2 * ζz^2) * ph1
    du[14] = (-0.5 / M) * k * f(t) * (du[10] - (K2 + safe_ddfoverf(t)) * ps0) + (Kp2 - ky^2) * php1 + (Km2 - ky^2) * phm1 - (0.5 * γ / M) * (ph2 - K2 * ph0)
    du[15] = (-0.5 * ky / M) * (f(t) * (m + 1) + 1) * (du[12] - (Kp2 + safe_ddfoverf(t)) * psp0) - (0.5 * γ / M) * (php2 - Kp2 * php0) + (K2 - ky^2 * ζz^2) * ph1
end

I12 = Matrix(I, 12, 12)

# coupling matrix (between ∂ₓ³ϕ's)
C = [0 1 0;
    1 0 1;
    0 1 0]

mass_matrix = zeros(15, 15)
mass_matrix[1:12, 1:12] = I12
mass_matrix[13:15, 13:15] = C

function bca!(res, u, p)
    # LEFT BOUNDARY (x=0)
    res[1] = u[10]      # ψm'(0)=0,   Dirichlet on right boundary (even)
    res[2] = u[7]       # ψm-1(0)=0, Dirichlet on right boundary (odd)
    res[3] = u[11]      # ψm+1(0)=0, Dirichlet on right boundary (odd)
    res[4] = u[3]       # ϕm(0)=0,    Dirichlet on left boundary (odd)
    res[5] = u[2]       # ϕm-1'(0)=0,  Dirichlet on left boundary (even)
    res[6] = u[6]       # ϕm+1'(0)=0,  Dirichlet on left boundary (even)
    res[7] = u[9] - 1   # ψm(0) = 1,  extra constraint to fix unknown parameter
end

function bcb!(res, u, p)
    # RIGHT BOUNDARY (x=L)
    res[1] = u[9]     # ψm(L)=0,    Dirichlet on right boundary
    res[2] = u[7]     # ψm-1(L)=0,  Dirichlet on right boundary
    res[3] = u[11]    # ψm+1(L)=0,  Dirichlet on right boundary
    res[4] = u[3]    # ϕm(L)=0,    Dirichlet on right boundary
    res[5] = u[1]    # ϕm-1(L)=0,  Dirichlet on right boundary
    res[6] = u[5]    # ϕm+1(L)=0,  Dirichlet on right boundary
    res[7] = u[4]    # ϕm'(L)=0,   Neumann on right boundary
    res[8] = u[2]    # ϕm-1'(L)=0, Neumann on right boundary
    res[9] = u[6]    # ϕm+1'(L)=0, Neumann on right boundary
end

function initial_guess(p, t)
    evenpsi = exp(-t^2)
    oddpsi = -2t * exp(-t^2)
    oddphi = t * exp(-t^2)
    evenphi = (1 - 2t^2) * exp(-t^2)

    return [
        evenphi;    # phm0, even
        oddphi;   # phm1, odd
        oddphi;    # ph0, odd
        evenphi;   # ph1, even
        evenphi;    # php0, even
        oddphi;   # php1, odd
        oddpsi;    # psm0, odd
        evenpsi;   # psm1, even
        evenpsi;    # ps0, even
        oddpsi;   # ps1, odd
        oddpsi;    # psp0, odd
        evenpsi;   # psp1, even
        evenphi;   # phm2, even
        oddphi;   # ph2, odd
        evenphi    # php2, even
    ]
end

# u0 = [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1];

fun = BVPFunction(tftearing!, (bca!, bcb!), mass_matrix=mass_matrix, twopoint=Val(true), bcresid_prototype=(zeros(7), zeros(9)))
prob = TwoPointBVProblem(fun, initial_guess, tspan, [γ_guess], fit_parameters=true)

sol = solve(prob, MIRK6(), dt=0.01,
    adaptive=true,
    progress=true,
    verbose=true,
    maxiters=500
)

# print the estimated value of Q which satisfies the BCs
println("γ fitted: ", sol.prob.p[1])

plot(sol, idxs=(0, 9), label=L"ψ(x)")
plot!(sol, idxs=(0, 3), label=L"φ(x)", xlabel="x", legend=:topright)

# plot(sol)