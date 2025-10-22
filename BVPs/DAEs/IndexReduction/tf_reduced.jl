using ModelingToolkit, OrdinaryDiffEq
using BoundaryValueDiffEq, LinearAlgebra, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel TWOFLUID begin
    @parameters begin
        # growth rate guess
        γ = 0.01

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
    end

    @variables begin
        ψm(t)
        ψ(t)
        ψp(t)
        φm(t)
        φ(t)
        φp(t)
        f(t)
    end

    @equations begin
        f ~ tanh(t)
        # m-1 mode equations
        γ * (D(D(φm)) - Km2 * φm) + (M / 2) * (D(D(D(φ))) - (Km2 - ζz^2 * ky^2) * D(φ)) ~ -ky * (f * (m - 1) - 1) * (D(D(ψm)) - Km2 * ψm - (D(D(f)) / f) * ψm)
        γ * ψm + (M / 2) * D(ψ) - ky * (f * (m - 1) - 1) * φm - (1 / S) * (D(D(ψm)) - Km2 * ψm) ~ 0

        # m mode equations
        γ * (D(D(φ)) - K2 * φ) + (M / 2) * (D(D(D(φp))) - (Kp2 - ky^2) * D(φp) + D(D(D(φm))) - (Kp2 - ky^2) * D(φm)) ~ -k * f * (D(D(ψ)) - K2 * ψ - (D(D(f)) / f) * ψ)
        γ * ψ + (M / 2) * (D(ψp) + D(ψm)) - k * f * φ - (1 / S) * (D(D(ψ)) - K2 * ψ) ~ 0

        # m+1 mode equations
        γ * (D(D(φp)) - Kp2 * φp) + (M / 2) * (D(D(D(φ))) - (Kp2 - ζz^2 * ky^2) * D(φ)) ~ -ky * (f * (m + 1) + 1) * (D(D(ψp)) - Kp2 * ψp - (D(D(f)) / f) * ψp)
        γ * ψp + (M / 2) * D(ψ) - ky * (f * (m + 1) + 1) * φp - (1 / S) * (D(D(ψp)) - Kp2 * ψp) ~ 0

    end
end

@mtkcompile tfmodel = TWOFLUID() # model building, index reduction, simplification in one step!

# now create the boundary value problem 
const L = 15.0
tspan = (0.0, L)

function bca!(res, u, p)
    # LEFT BOUNDARY (x=0)
    # 1) ψm'(0) = 0  -> ψm_t is index 7
    res[1] = u[7]

    # 2) ψ_{m-1}(0) = 0 -> ψm (psi minus) is index 5
    res[2] = u[5]

    # 3) ψ_{m+1}(0) = 0 -> ψp is index 1
    res[3] = u[1]

    # 4) φ_{m}(0) = 0 -> φ (phi center) is index 9
    res[4] = u[9]

    # 5) φ_{m-1}'(0) = 0 -> φm_t is index 13
    res[5] = u[13]

    # 6) φ_{m+1}'(0) = 0 -> φp_t is index 12
    res[6] = u[12]

    # 7) ψ_{m}(0) = 1  -> ψ (psi center) is index 3
    res[7] = u[3] - 1
end

function bcb!(res, u, p)
    # RIGHT BOUNDARY (x=L)
    # 1) ψ_{m}(L) = 0 -> ψ is index 3
    res[1] = u[3]

    # 2) ψ_{m-1}(L) = 0 -> ψm is index 5
    res[2] = u[5]

    # 3) ψ_{m+1}(L) = 0 -> ψp is index 1
    res[3] = u[1]

    # 4) φ_{m}(L) = 0 -> φ is index 9
    res[4] = u[9]

    # 5) φ_{m-1}(L) = 0 -> φm is index 8
    res[5] = u[8]

    # 6) φ_{m+1}(L) = 0 -> φp is index 2
    res[6] = u[2]

    # 7) φ_{m}'(L) = 0 -> φ_t is index 10
    res[7] = u[10]

    # 8) φ_{m-1}'(L) = 0 -> φm_t is index 13
    res[8] = u[13]

    # 9) φ_{m+1}'(L) = 0 -> φp_t is index 12
    res[9] = u[12]
end

function initial_guess(p, t)
    # Base building blocks
    evenpsi = exp(-t^2)
    oddpsi = -2t * exp(-t^2)
    evenphi = (1 - 2t^2) * exp(-t^2)
    oddphi = t * exp(-t^2)

    # Their simple “derivative-like” versions
    evenpsi_t = oddpsi          # derivative of evenψ is odd
    oddpsi_t = 2 * (t^2 - 1) * exp(-t^2)   # derivative of oddψ is even
    evenphi_t = -2t * (3 - 2t^2) * exp(-t^2)  # derivative of evenφ is odd
    oddphi_t = (1 - 2t^2) * exp(-t^2)         # derivative of oddφ is even

    # For higher derivatives, we’ll keep same parity trend but simplified (still smooth)
    evenphi_tt = evenphi   # even stays even for 2nd derivative
    evenphi_ttt = oddphi    # odd for 3rd derivative
    oddphi_tt = evenphi
    oddphi_ttt = oddphi

    # Now assign according to your new state vector
    return [
        evenpsi;     # 1: ψp(t)       (even)
        evenphi;     # 2: φp(t)       (even)
        evenpsi;     # 3: ψ(t)        (even)
        oddpsi;      # 4: ψp_t(t)     (odd)
        evenpsi;     # 5: ψm(t)       (even)
        oddpsi;      # 6: ψ_t(t)      (odd)
        oddpsi;      # 7: ψm_t(t)     (odd)
        evenphi;     # 8: φm(t)       (even)
        evenphi;     # 9: φ(t)        (even)
        oddphi;      #10: φ_t(t)      (odd)
        evenphi_tt;  #11: φ_tt(t)     (even)
        oddphi_t;    #12: φp_t(t)     (odd)
        oddphi_t;    #13: φm_t(t)     (odd)
        evenphi_ttt; #14: φ_ttt(t)    (odd)
        oddphi_ttt   #15: φm_ttt(t)   (even-ish/placeholder)
    ]
end

# Convert MTK model to numerical residual
odeprob = ODEProblem(tfmodel, [], tspan, guesses = [γ = 0.01])
# mass_matrix = calculate_massmatrix(tfmodel)

# Construct Boundary Value Problem


#=
prob = ODEroblem(tfmodel,
    bc=(bca!, bcb!),
    initial_guess, tspan, [γ_guess],
    mass_matrix=mass_matrix,
    twopoint=Val(true),
    bcresid_prototype=(zeros(7), zeros(9), fit_parameters=true)
)

sol = solve(prob, GeneralMIRK4(), dt=0.1, maxiters=100, verbose=true)

plot(sol)
println("Estimated γ = ", sol.param_estimates[1])
=#
