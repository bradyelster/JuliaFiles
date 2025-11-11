using Plots, LaTeXStrings

L = 15.0
f(t) = tanh(t)

# mode numbers of interest
m = 2
n = 1

mp1 = m + 1
mm1 = m - 1

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

Bz = k / kn # chosen such that kₘ ⋅ Β₀ = 0 at x=0

# Define domain
tvals = LinRange(-L, L, 1000)

# Define the functions to plot
f_mplus(t) = kp * (1 + 0.1 * f(t)) - kn * Bz
f_mminus(t) = km * (1 + 0.1 * f(t)) - kn * Bz
f_m(t) = k * (1 + 0.1 * f(t)) - kn * Bz

# Plot the functions
plot(tvals, f_mplus.(tvals), label="$mp1/$n", linewidth=2)
plot!(tvals, f_m.(tvals), label="$m/$n (res.)", linewidth=2)
plot!(tvals, f_mminus.(tvals), label="$mm1/$n", linewidth=2)

xlabel!("x")
# ylabel!("f(t)")
title!("Resonance Conditions")