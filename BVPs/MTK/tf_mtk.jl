using ModelingToolkit, BoundaryValueDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

@constants begin
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
    M = 1 / 2 #S^(-1 / 2)
end

@parameters begin
    γ
end

@variables begin
    ψm(..)
    ψ(..)
    ψp(..)
    φm(..)
    φ(..)
    φp(..)
    f(..)
end;

eqns = [
    f(t) ~ tanh(t),
    # m-1 mode equations
    γ * (D(D(φm(t))) - Km2 * φm(t)) + (M / 2) * (D(D(D(φ(t)))) - (Km2 - ζz^2 * ky^2) * D(φ(t))) ~ -ky * ((f(t) * (m - 1) - 1) * D(D(ψm(t))) - Km2 * (f(t) * (m - 1) - 1) * ψm(t) - (D(D(f(t))) * ((m - 1) - 1)) * ψm(t)),
    γ * ψm(t) + (M / 2) * D(ψ(t)) - ky * (f(t) * (m - 1) - 1) * φm(t) ~ (1 / S) * (D(D(ψm(t))) - Km2 * ψm(t)),

    # m mode equations
    γ * (D(D(φ(t))) - K2 * φ(t)) + (M / 2) * (D(D(D(φp(t)))) - (Kp2 - ky^2) * D(φp(t)) + D(D(D(φm(t)))) - (Km2 - ky^2) * D(φm(t))) ~ -k * f(t) * D(D(ψ(t))) + k * K2 * f(t) * ψ(t) + k * D(D(f(t))) * ψ(t),
    γ * ψ(t) + (M / 2) * (D(ψp(t)) + D(ψm(t))) - k * f(t) * φ(t) - (1 / S) * (D(D(ψ(t))) - K2 * ψ(t)) ~ 0,

    # m+1 mode equations
    γ * (D(D(φp(t))) - Kp2 * φp(t)) + (M / 2) * (D(D(D(φ(t)))) - (Kp2 - ζz^2 * ky^2) * D(φ(t))) ~ -ky * ((f(t) * (m + 1) + 1) * D(D(ψp(t))) - Kp2 * (f(t) * (m + 1) + 1) * ψp(t) - (D(D(f(t))) * ((m + 1) + 1)) * ψp(t)),
    γ * ψp(t) + (M / 2) * D(ψ(t)) - ky * (f(t) * (m + 1) + 1) * φp(t) ~ (1 / S) * (D(D(ψp(t))) - Kp2 * ψp(t))
];

global L = 15.0

# right boundary conditions
constr = [
    # all functions zero at right side (6)
    ψm(L) ~ 0.0,
    φm(L) ~ 0.0,
    ψ(L) ~ 0.0,
    φ(L) ~ 0.0,
    ψp(L) ~ 0.0,
    φp(L) ~ 0.0,
];

@mtkcompile bvp = System(
    eqns, t;
    constraints=constr
);

prob = BVProblem(
    bvp,
    [
        # all functions zero at right side (6 + 6 = 12)
        ψm(t) => 0.0,
        φm(t) => 0.0,
        ψ(t) => 0.0,
        φ(t) => 0.0,
        ψp(t) => 0.0,
        φp(t) => 0.0,

        # try 2 more to get to 14 BCs (12 + 2 = 14) ✅
        D(φ(t)) => 0.0,
        D(ψ(t)) => 0.0
    ],
    (-L, L),
    guesses=[
        γ => 0.1,
        D(ψp(t)) => 0.1,
        D(ψm(t)) => 0.1,
        D(D(φ(t))) => 0.1,
        D(φp(t)) => 0.1,
        D(φm(t)) => 0.1,
        D(D(D(φ(t)))) => 0.1,
        D(D(D(φm(t)))) => 0.1
    ],
    #jac=true,
    #sparse=true
)

solve(prob, MIRK4(), dt=0.05)

#plot(sol)