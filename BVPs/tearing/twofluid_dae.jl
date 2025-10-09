# Tearing Mode - Numerical Solution, 10/03/2025
# Two-Fluid Tearing Equations, DAE Solution

using BoundaryValueDiffEq, LinearAlgebra, Plots, LaTeXStrings
L = 12.0
tspan = (0.0, L)
f(t) = tanh(t)
ddf(t) = -2 * tanh(t) * sech(t)^2

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
Km2 = km^2 + kn^2
K2 = k^2 + kn^2
Kp2 = kp^2 + kn^2

# Zeta-components (CHANGE THESE LATER)
ζz = 1
ζy = 1

# Lundquist number
S = 100
ε = S^(-2/3)

# Growth rate
Q_guess = 0.267

# Alfvénic Mach number
M = S^(-1/2)

function tftearing!(du, u, p, t)
    phm0, phm1, 
    ph0, ph1, 
    php0, php1, 
    psm0, psm1, 
    ps0, ps1, 
    psp0, psp1, 
    phm2, ph2, php2 = u
    Q = p[1]
    du[1] = phm1
    du[2] = phm2
    du[3] = ph1
    du[4] = ph2
    du[5] = php1
    du[6] = php2
    du[7] = psm1
    du[8] = (Q/ε + Km2)*psm0 + (M/(2*ε^(3/2)))*ps1 - (ky*(f(t)*(m-1)-1)/ε)*phm0 # ψm-1''
    du[9] = ps1
    du[10] = (Q/ε + K2)*ps0 + (M/(2*ε^(3/2)))*(psp1 + psm1) - (k*f(t)/ε)*ph0 # ψm''
    du[11] = psp1
    du[12] = (Q/ε + Kp2)*psp0 + (M/(2*ε^(3/2)))*ps1 - (ky*(f(t)*(m+1)+1)/ε)*php0 # ψm+1''
    du[13] = -(2/(M*ε^(1/2)))*ky*(f(t)*(m-1)-1)*(du[8]-Km2*psm0-(ddf(t)/f(t))*psm0) - ((2*ε^(1/2))/M)*Q*(phm2 - Km2*phm0) + (K2 - ζz^2*ky^2)*ph1
    du[14] = -(2/(M*ε^(1/2)))*k*f(t)*(du[10]-K2*ps0-(ddf(t)/f(t))*ps0) - ((2*ε^(1/2))/M)*Q*(ph2 - K2*ph0) + (Kp2 - ky^2)*php1 + (Km2 - ky^2)*phm1
    du[15] = -(2/(M*ε^(1/2)))*ky*(f(t)*(m+1)+1)*(du[12]-Kp2*psp0-(ddf(t)/f(t))*psp0) - ((2*ε^(1/2))/M)*Q*(php2 - Kp2*php0) + (K2 - ζz^2*ky^2)*ph1
end

I12 = Matrix(I, 12, 12)

# coupling matrix (between ∂ₓ³ϕ's)
C = [0 1 0;
     1 0 1;
     0 1 0]

mass_matrix = zeros(15, 15)
# place the identity block

mass_matrix[1:12, 1:12] = I12

# place the singular coupling matrix block
mass_matrix[13:15, 13:15] = C

function bc!(res, u, p, t)
    Q = p[1]
    # note: u[time][variable index] is the format (I think!)
    # LEFT BOUNDARY (x=0)
    res[1] = u[1][10]      # ψm'(0)=0,   Dirichlet on right boundary (even)
    res[2] = u[1][8]       # ψm-1'(0)=0, Dirichlet on right boundary (even)
    res[3] = u[1][12]      # ψm+1'(0)=0, Dirichlet on right boundary (even)
    res[4] = u[1][3]       # ϕm(0)=0,    Dirichlet on left boundary (odd)
    res[5] = u[1][1]       # ϕm-1(0)=0,  Dirichlet on left boundary (odd)
    res[6] = u[1][5]       # ϕm+1(0)=0,  Dirichlet on left boundary (odd)
    res[7] = u[1][9] - 1  # ψm(0) = 1:  extra constraint to fix unknown parameter Q
    
    # RIGHT BOUNDARY (x=L)
    res[8] = u[end][9]     # ψm(L)=0,    Dirichlet on right boundary
    res[9] = u[end][7]     # ψm-1(L)=0,  Dirichlet on right boundary
    res[10] = u[end][11]    # ψm+1(L)=0,  Dirichlet on right boundary
    res[11] = u[end][3]    # ϕm(L)=0,    Dirichlet on right boundary
    res[12] = u[end][1]    # ϕm-1(L)=0,  Dirichlet on right boundary
    res[13] = u[end][5]    # ϕm+1(L)=0,  Dirichlet on right boundary
    res[14] = u[end][4]    # ϕm'(L)=0,   Neumann on right boundary
    res[15] = u[end][2]    # ϕm-1'(L)=0, Neumann on right boundary
    res[16] = u[end][6]    # ϕm+1'(L)=0, Neumann on right boundary
end

function initial_guess(p, t)
    ψ  = exp(-t^2)                   
    ψ1 = -2t * exp(-t^2)             
    ϕ  = t * exp(-t^2)               
    ϕ1 = (1 - 2t^2) * exp(-t^2)
    ϕ2 = t * (-6 + 4t^2) * exp(-t^2)

    return [
        ϕ;    # phm0
        ϕ1;   # phm1
        ϕ;    # ph0
        ϕ1;   # ph1
        ϕ;    # php0
        ϕ1;   # php1
        ψ;    # psm0
        ψ1;   # psm1
        ψ;    # ps0
        ψ1;   # ps1
        ψ;    # psp0
        ψ1;   # psp1
        ϕ2;   # phm2
        ϕ2;   # ph2
        ϕ2    # php2
    ]
end

fun = BVPFunction(tftearing!, bc!, mass_matrix = mass_matrix, bcresid_prototype=zeros(16))
prob = BVProblem(fun, initial_guess, tspan, [Q_guess], fit_parameters=true)
sol = solve(prob, MIRK4(), dt=0.01, save_everystep=false, verbose=true, maxiters=100)

# print the estimated value of Q which satisfies the BCs
println(sol.prob.p[1])

plot(sol, idxs=(0, 9), label=L"ψ(x)", continuity=:right)
plot!(sol, idxs=(0, 3), label=L"φ(x)", xlabel="x", legend=:topright, continuity=:right)