using ApproxFun
using LinearAlgebra
using Plots

# ---------------------------
# Problem parameters
# ---------------------------
L = 20.0
ε = 0.5      # epsilon (same as k before)
η = 0.001    # eta = 1/S, so if S=1000, η=0.001

f(x) = tanh(x)
df(x) = sech(x)^2
ddf(x) = -2 * tanh(x) * sech(x)^2

# ---------------------------
# Function space
# ---------------------------
C = Chebyshev(0 .. L)
𝒟 = Derivative(C)
B = Dirichlet(C)

# Create Fun objects
x = Fun(identity, C)
f_fun = f(x)
df_fun = df(x)
ddf_fun = ddf(x)

# ---------------------------
# Build generalized eigenvalue problem: L*u = λ*M*u
# where u = [ψ; ϕ]
# ---------------------------

# For single field with Dirichlet BCs, following the quantum example:
# We need to reformulate as a single operator problem or use a different approach

# Let's try converting to standard eigenvalue by "inverting" M
# M*u_λ = λ*L*u → solve M\L*u = (1/λ)*u

# LEFT-HAND SIDE (L operator) - with BCs incorporated
Lψψ = η * (𝒟^2 - ε^2 * I)
Lψϕ = -f_fun * I
Lϕψ = -f_fun * (𝒟^2 - ε^2 * I) + ddf_fun * I
Lϕϕ = 0 * I

# RIGHT-HAND SIDE (M operator)
Mψψ = I
Mψϕ = 0 * I
Mϕψ = 0 * I
Mϕϕ = 𝒟^2 - ε^2 * I

# Create the combined operator for standard eigenvalue problem
# We'll solve: (M^-1 * L) * u = λ * u
# This means we need L_effective = M\L

# Build block systems
Lblock = [Lψψ Lψϕ;
    Lϕψ Lϕϕ]

Mblock = [Mψψ Mψϕ;
    Mϕψ Mϕϕ]

# ---------------------------
# Solve as standard eigenvalue problem
# ---------------------------
nev = 20  # number of eigenvalues to compute

println("Solving eigenvalue problem with:")
println("  Domain: [0, $L]")
println("  ε = $ε")
println("  η = $η (S = $(1/η))")
println("  Boundary conditions: ψ(0)=ψ(L)=0, ϕ(0)=ϕ(L)=0")
println()

try
    # Solve M\L as standard eigenvalue problem with BCs
    # The eigenvalues will be λ from the original generalized problem
    L_effective = Mblock \ Lblock

    λ, v = ApproxFun.eigs(B, L_effective, nev, tolerance=1E-8)

    # Find index of largest real eigenvalue
    max_idx = argmax(real.(λ))
    λ_max = λ[max_idx]

    println("Largest real eigenvalue:")
    println("λ_max = $(round(real(λ_max), digits=6)) + $(round(imag(λ_max), digits=6))im")
    println()

    # Extract corresponding eigenfunctions
    ψ_max = v[max_idx][1]
    ϕ_max = v[max_idx][2]

    # ---------------------------
    # Plot eigenfunctions
    # ---------------------------
    xvals = range(0, L, length=500)

    # Evaluate at points
    ψvals = ψ_max.(xvals)
    ϕvals = ϕ_max.(xvals)

    # Create plot
    p = plot(xlabel="x", ylabel="Amplitude",
        title="Eigenfunctions for λ = $(round(real(λ_max), digits=6))",
        legend=:best, size=(900, 600))

    plot!(p, xvals, real.(ψvals), label="ψ(x)", linewidth=2, color=:blue)
    plot!(p, xvals, real.(ϕvals), label="ϕ(x)", linewidth=2, color=:red, linestyle=:dash)
    plot!(p, [0, L], [0, 0], color=:black, linestyle=:dot, label="", alpha=0.3)

    display(p)

    # Verify boundary conditions
    println("Boundary condition verification:")
    println("ψ(0) = $(ψ_max(0.0)) (should be ≈ 0)")
    println("ψ(L) = $(ψ_max(L)) (should be ≈ 0)")
    println("ϕ(0) = $(ϕ_max(0.0)) (should be ≈ 0)")
    println("ϕ(L) = $(ϕ_max(L)) (should be ≈ 0)")

catch e
    println("Error with eigs: $e")
    println("\nStacktrace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end