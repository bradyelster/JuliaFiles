using BoundaryValueDiffEq, Plots, LaTeXStrings

L = 15.0
tspan = (-L, L)
F(t) = tanh(t)
Fpp(t) = -2 * tanh(t) * sech(t)^2

k = 0.5

function tearing!(du, u, p, t)
    ψ, ψp, φ, φp = u
    S, γ = p
    du[1] = ψp
    du[2] = k^2 * ψ + S * γ * (ψ - F(t) * φ)
    du[3] = φp
    du[4] = (k^2 + S * F(t)^2 / γ) * φ + (Fpp(t) / γ^2 - S * F(t) / γ) * ψ
end

function bca!(res, u, p)
    res[1] = u[1]  # ψ(-L)=0
    res[2] = u[3]  # φ(-L)=0
end

function bcb!(res, u, p)
    res[1] = u[1]  # ψ(L)=0
    res[2] = u[3]  # φ(L)=0
end

guess(p, t) = [exp(-t^2); -2 * t * exp(-t^2); t * exp(-t^2); (1 - 2 * t^2) * exp(-t^2)]

# Continuation method: Start with small S, gradually increase to target S=100
println("Starting continuation method...")
S_values = [1.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0]
γ_initial = 0.16

sol = guess  # Start with initial guess function

for (i, S_target) in enumerate(S_values)
    println("\nIteration $i: S = $S_target")
    
    # Create BVP with current S value, fitting γ
    bvp = TwoPointBVProblem(tearing!, (bca!, bcb!), sol, tspan, [S_target, γ_initial],
        bcresid_prototype=(zeros(2), zeros(2)), fit_parameters=[false, true])
    
    # Solve with adaptive tolerances (tighter as S increases)
    abstol = max(1e-8, 1e-6 / sqrt(S_target))
    reltol = max(1e-6, 1e-4 / sqrt(S_target))
    
    global sol = solve(bvp, MIRK6(), dt=0.05, abstol=abstol, reltol=reltol, maxiters=1000)
    
    # Update γ guess for next iteration
    global γ_initial = sol.prob.p[2]
    println("  Fitted γ = $γ_initial")
    
    # Use the solution as the next guess (convert to function form)
    global sol = (p, t) -> sol(t)
end

println("\n" * "="^50)
println("Final solution at S = 100:")
println("Fitted γ = ", sol.u.p[2])
println("="^50)

# Plot the final solution
p1 = plot(sol, idxs=1, xlabel=L"t", ylabel=L"\psi(t)", title="Wave function ψ (S=100)", lw=2, legend=false)
p2 = plot(sol, idxs=2, xlabel=L"t", ylabel=L"\psi'(t)", title="Derivative ψ' (S=100)", lw=2, legend=false)
p3 = plot(sol, idxs=3, xlabel=L"t", ylabel=L"\phi(t)", title="Wave function φ (S=100)", lw=2, legend=false)
p4 = plot(sol, idxs=4, xlabel=L"t", ylabel=L"\phi'(t)", title="Derivative φ' (S=100)", lw=2, legend=false)
plot(p1, p2, p3, p4, layout=(2,2), size=(800,600))

# Optional: Plot evolution of γ with S
# S_hist = [1.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0]
# γ_hist = [...] # Store γ values from each iteration
# plot(S_hist, γ_hist, xlabel="S", ylabel="γ", marker=:o, 
#      title="Growth rate vs Lundquist number", xscale=:log10, legend=false)