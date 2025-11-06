using Plots, LaTeXStrings

L = 15.0
f(t) = tanh(t)
# g(t) = -2 * sech(t)^2

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

Bz = k / kn

# Lundquist number
S = 10

# Alfvénic Mach number
M = 1 / 2 #S^(-1 / 2)

# Define domain
tvals = LinRange(-L, L, 1000)

# Define the functions to plot
f_mplus(t) = kp * (1 + f(t)) - kn * Bz
f_mminus(t) = km * (1 + f(t)) - kn * Bz
f_m(t) = k * (1 + f(t)) - kn * Bz

# Plot the functions
plot(tvals, f_mplus.(tvals), label="m+1")
plot!(tvals, f_m.(tvals), label="m")
plot!(tvals, f_mminus.(tvals), label="m-1")

# Compute and plot the values at t = 0
t0 = 0.0
points_y = [f_mplus(t0), f_mminus(t0), f_m(t0)]

scatter!([t0, t0, t0], points_y, markershape=:circle, label=false)

xlabel!("t")

# measure vertical distance between points 
Δ_plus = f_mplus(t0) - f_m(t0)
Δ_minus = f_m(t0) - f_mminus(t0)

println("Δ_plus  = ", Δ_plus)
println("Δ_minus = ", Δ_minus)

println("ky = ", ky)
println("Δ_plus ≈ ky ? ", isapprox(Δ_plus, ky; atol=1e-8))
println("Δ_minus ≈ ky ? ", isapprox(Δ_minus, ky; atol=1e-8))
