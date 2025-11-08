using Plots, LaTeXStrings

L = 15.0
f(t) = tanh(t)

# mode numbers of interest
m = 3
n = 2

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

Bz = k / kn

# Lundquist number
S = 10

# Alfvénic Mach number
M = 1 / 2

# Define domain
tvals = LinRange(-L, L, 1000)

# Define the functions to plot
f_mplus(t) = kp * (1 + f(t)) - kn * Bz
f_mminus(t) = km * (1 + f(t)) - kn * Bz
f_m(t) = k * (1 + f(t)) - kn * Bz

# Bisection method to find zeros
function find_zero(func, a, b; tol=1e-10, max_iter=100)
    fa = func(a)
    fb = func(b)
    
    # Check if there's a sign change
    if fa * fb > 0
        return nothing
    end
    
    for i in 1:max_iter
        c = (a + b) / 2
        fc = func(c)
        
        if abs(fc) < tol || (b - a) / 2 < tol
            return c
        end
        
        if fa * fc < 0
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end
    end
    
    return (a + b) / 2
end

# Find zeros for each function
zeros_mplus = []
zeros_m = []
zeros_mminus = []

# Search for zeros in overlapping intervals
interval_size = 1.0
for t_start in -L:interval_size:(L-interval_size)
    t_end = t_start + interval_size
    
    # Check each function
    z = find_zero(f_mplus, t_start, t_end)
    if z !== nothing
        push!(zeros_mplus, z)
    end
    
    z = find_zero(f_m, t_start, t_end)
    if z !== nothing
        push!(zeros_m, z)
    end
    
    z = find_zero(f_mminus, t_start, t_end)
    if z !== nothing
        push!(zeros_mminus, z)
    end
end

# Plot the functions
plot(tvals, f_mplus.(tvals), label="m+1", linewidth=2)
plot!(tvals, f_m.(tvals), label="m", linewidth=2)
plot!(tvals, f_mminus.(tvals), label="m-1", linewidth=2)

# Plot the zeros
scatter!(zeros_mplus, zeros(length(zeros_mplus)), 
         markershape=:circle, markersize=6, label="zeros m+1", color=1)
scatter!(zeros_m, zeros(length(zeros_m)), 
         markershape=:circle, markersize=6, label="zeros m", color=2)
scatter!(zeros_mminus, zeros(length(zeros_mminus)), 
         markershape=:circle, markersize=6, label="zeros m-1", color=3)

xlabel!("t")
ylabel!("f(t)")
title!("Functions and their zeros")

# Print the zeros
println("Zeros of f_mplus:  ", zeros_mplus)
println("Zeros of f_m:      ", zeros_m)
println("Zeros of f_mminus: ", zeros_mminus)
