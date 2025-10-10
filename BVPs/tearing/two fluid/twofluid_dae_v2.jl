# Tearing Mode - Numerical Solution, 10/03/2025
# Two-Fluid Tearing Equations, DAE Solution - OPTIMIZED

using BoundaryValueDiffEq, LinearAlgebra, Plots, LaTeXStrings
using StaticArrays  # Add for performance

L = 12.0
tspan = (0.0, L)

# Use @inline for small functions
@inline f(t) = tanh(t)
@inline ddf(t) = -2 * tanh(t) * sech(t)^2

# mode numbers of interest
m, n = 2, 1

# device parameters
R0 = 1

# poloidal wave numbers 
ky = 0.5
km = (m-1)*ky
k = (m)*ky
kp = (m+1)*ky

# toroidal wave number
kz = 1/R0
kn = (n)*kz

# combined wave numbers
const Km2 = km^2 + kn^2
const K2 = k^2 + kn^2
const Kp2 = kp^2 + kn^2

# Zeta-components (CHANGE THESE LATER)
ζz = 1
ζy = 1

# Lundquist number
S = 100
ε = S^(-2/3)

# Pre-compute constants
const inv_ε = 1/ε
const inv_ε2 = 1/ε^2
const inv_ε_sqrt = 1/sqrt(ε)
const ε_sqrt = sqrt(ε)

# Growth rate
Q_guess = 0.267

# Alfvénic Mach number
M = S^(-1/2)

# Pre-compute more constants
const M_inv_ε32 = M/(2*ε^(3/2))
const two_inv_M_ε_sqrt = 2/(M*sqrt(ε))
const two_ε_sqrt_inv_M = (2*sqrt(ε))/M
const ky_m_minus_1 = ky*(m-1)
const ky_m_plus_1 = ky*(m+1)
const K2_minus_ζz2_ky2 = K2 - ζz^2*ky^2
const Kp2_minus_ky2 = Kp2 - ky^2
const Km2_minus_ky2 = Km2 - ky^2

function tftearing!(du, u, p, t)
    phm0, phm1, 
    ph0, ph1, 
    php0, php1, 
    psm0, psm1, 
    ps0, ps1, 
    psp0, psp1, 
    phm2, ph2, php2 = u
    Q = p[1]
    
    # Cache f(t) and ddf(t) to avoid recomputation
    ft = f(t)
    ddft = ddf(t)
    
    # Pre-compute common terms
    Q_inv_ε = Q*inv_ε
    ddft_over_ft = ddft/ft
    ky_ft_m_minus_1_minus_1 = ky_m_minus_1*ft - ky
    ky_ft_m_plus_1_plus_1 = ky_m_plus_1*ft + ky
    k_ft = k*ft
    
    du[1] = phm1
    du[2] = phm2
    du[3] = ph1
    du[4] = ph2
    du[5] = php1
    du[6] = php2
    du[7] = psm1
    du[8] = (Q_inv_ε + Km2)*psm0 + M_inv_ε32*ps1 - (ky_ft_m_minus_1_minus_1*inv_ε)*phm0 # ψm-1''
    du[9] = ps1
    du[10] = (Q_inv_ε + K2)*ps0 + M_inv_ε32*(psp1 + psm1) - k_ft*inv_ε*ph0 # ψm''
    du[11] = psp1
    du[12] = (Q_inv_ε + Kp2)*psp0 + M_inv_ε32*ps1 - (ky_ft_m_plus_1_plus_1*inv_ε)*php0 # ψm+1''
    
    # Cache du[8], du[10], du[12] terms for reuse
    du8_term = du[8] - Km2*psm0 - ddft_over_ft*psm0
    du10_term = du[10] - K2*ps0 - ddft_over_ft*ps0
    du12_term = du[12] - Kp2*psp0 - ddft_over_ft*psp0
    
    du[13] = -two_inv_M_ε_sqrt*ky_ft_m_minus_1_minus_1*du8_term - two_ε_sqrt_inv_M*Q*(phm2 - Km2*phm0) + K2_minus_ζz2_ky2*ph1
    du[14] = -two_inv_M_ε_sqrt*k_ft*du10_term - two_ε_sqrt_inv_M*Q*(ph2 - K2*ph0) + Kp2_minus_ky2*php1 + Km2_minus_ky2*phm1
    du[15] = -two_inv_M_ε_sqrt*ky_ft_m_plus_1_plus_1*du12_term - two_ε_sqrt_inv_M*Q*(php2 - Kp2*php0) + K2_minus_ζz2_ky2*ph1
end

const I12 = Matrix(I, 12, 12)

# coupling matrix (between ∂ₓ³ϕ's)
const C = [0 1 0;
           1 0 1;
           0 1 0]

const mass_matrix = let mm = zeros(15, 15)
    # place the identity block
    mm[1:12, 1:12] = I12
    # place the singular coupling matrix block
    mm[13:15, 13:15] = C
    mm
end

function bc!(res, u, p, t)
    Q = p[1]
    # note: u[time][variable index] is the format (I think!)
    # LEFT BOUNDARY (x=0)
    res[1] = u[1][10]      # ψm'(0)=0,   Dirichlet on right boundary (even)
    res[2] = u[1][8]       # ψm-1'(0)=0, Dirichlet on right boundary (even)
    res[3] = u[1][12]      # ψm+1'(0)=0, Dirichlet on right boundary (even)
    res[4] = u[1][3]       # ϕm(0)=0,    Dirichlet on left boundary (odd)
    #res[5] = u[1][1]       # ϕm-1(0)=0,  Dirichlet on left boundary (odd)
    #res[6] = u[1][5]       # ϕm+1(0)=0,  Dirichlet on left boundary (odd)
    res[5] = u[1][2]       # ϕm-1'(0)=0,  Dirichlet on left boundary (odd)
    res[6] = u[1][6]       # ϕm+1'(0)=0,  Dirichlet on left boundary (odd)
    res[7] = u[1][9] - 1   # ψm(0) = 1:  extra constraint to fix unknown parameter Q
    
    # RIGHT BOUNDARY (x=L)
    res[8] = u[end][9]     # ψm(L)=0
    res[9] = u[end][7]     # ψm-1(L)=0
    res[10] = u[end][11]   # ψm+1(L)=0
    res[11] = u[end][3]    # ϕm(L)=0
    res[12] = u[end][1]    # ϕm-1(L)=0
    res[13] = u[end][5]    # ϕm+1(L)=0
    res[14] = u[end][4]    # ϕm'(L)=0
    res[15] = u[end][2]    # ϕm-1'(L)=0
    res[16] = u[end][6]    # ϕm+1'(L)=0
end

function initial_guess(p, t)
    exp_t2 = exp(-t^2)
    t_exp = t * exp_t2
    t2 = t^2
    
    ψ  = exp_t2                   
    ψ1 = -2t * exp_t2             
    ϕ  = t_exp               
    ϕ1 = (1 - 2t2) * exp_t2
    ϕ2 = t * (-6 + 4t2) * exp_t2

    return SVector(
        ϕ,    # phm0
        ϕ1,   # phm1
        ϕ,    # ph0
        ϕ1,   # ph1
        ϕ,    # php0
        ϕ1,   # php1
        ψ,    # psm0
        ψ1,   # psm1
        ψ,    # ps0
        ψ1,   # ps1
        ψ,    # psp0
        ψ1,   # psp1
        ϕ2,   # phm2
        ϕ2,   # ph2
        ϕ2    # php2
    )
end

fun = BVPFunction(tftearing!, bc!, mass_matrix = mass_matrix, bcresid_prototype=zeros(16))
prob = BVProblem(fun, initial_guess, tspan, [Q_guess], fit_parameters=true)

# Try coarser initial dt
@time sol = solve(prob, MIRK4(), dt=0.01, abstol=1e-6, reltol=1e-4)

# print the estimated value of Q which satisfies the BCs
#println(sol.prob.p[1])

#plot(sol, idxs=(0, 9), label=L"ψ(x)", continuity=:right)
#plot!(sol, idxs=(0, 3), label=L"φ(x)", xlabel="x", legend=:topright, continuity=:right)