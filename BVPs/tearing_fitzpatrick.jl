using BoundaryValueDiffEq, Plots, LaTeXStrings

L = 15.0
tspan = (-L, L)
F(x) = tanh(x)
Fpp(x) = -2*tanh(x)*sech(x)^2

function system!(du,u,p,x)
    ψ, ψp, φ, φp = u
    S, k, γ = p
    du[1] = ψp
    du[2] = k^2*ψ + S*γ*(ψ - F(x)*φ)
    du[3] = φp
    du[4] = (k^2 + S*F(x)^2/γ)*φ + (Fpp(x)/γ^2 - S*F(x)/γ)*ψ
end

function bca!(res,u,p)
    res[1] = u[1]  # ψ(-L)=0
    res[2] = u[3]  # φ(-L)=0
end

function bcb!(res,u,p)
    res[1] = u[1]  # ψ(L)=0
    res[2] = u[3]  # φ(L)=0
end

function guess(p)
    # Return a function that takes x and returns the initial guess
    function initial_guess(x)
        S, k, γ = p
        
        # Handle both scalar and array inputs
        if x isa Number
            xi = x
            
            # ψ is even: ψ(-x) = ψ(x)
            # φ is odd: φ(-x) = -φ(x)
            
            decay = exp(-k*abs(xi))
            signx = xi == 0 ? 0.0 : sign(xi)
            
            # Even function for ψ
            ψ = decay * sech(xi)^2
            # Derivative of even function is odd
            ψp = -decay * 2*sech(xi)^2*tanh(xi) - k*signx*ψ
            
            # Odd function for φ: multiply by sign(x) or use tanh(x)
            φ = decay * tanh(xi) * sech(xi)  # tanh is odd, sech is even
            # Derivative of odd function is even
            φp = decay * (sech(xi)^3 - tanh(xi)^2*sech(xi)) - k*signx*φ
            
            return [ψ, ψp, φ, φp]
        else
            # x is an array
            u0 = zeros(4, length(x))
            
            for (i, xi) in enumerate(x)
                decay = exp(-k*abs(xi))
                signx = xi == 0 ? 0.0 : sign(xi)
                
                # Even function for ψ
                ψ = decay * sech(xi)^2
                ψp = -decay * 2*sech(xi)^2*tanh(xi) - k*signx*ψ
                
                # Odd function for φ
                φ = decay * tanh(xi) * sech(xi)
                φp = decay * (sech(xi)^3 - tanh(xi)^2*sech(xi)) - k*signx*φ
                
                u0[1, i] = ψ
                u0[2, i] = ψp
                u0[3, i] = φ
                u0[4, i] = φp
            end
            
            return u0
        end
    end
    
    return initial_guess
end

# Define parameters (you need to set these values)
S = 100      # Set your value for S
k = 0.5      # Set your value for k
γ_init = S^(-2/5) # Set your initial guess for γ

p = [S, k, γ_init]

function bca!(res,u,p)
    res[1] = u[1]  # ψ(-L)=0
    res[2] = u[3]  # φ(-L)=0
end

function bcb!(res,u,p)
    res[1] = u[1]   # ψ(L)=0
    res[2] = u[3]   # φ(L)=0
    res[3] = u[4]   # φ'(L)=0 (odd parity: φ'(0)=finite, but φ'(±L)≈0 for large L)
end

bvp = TwoPointBVProblem(system!, (bca!, bcb!), guess(p), tspan, p, 
                        bcresid_prototype=(zeros(2), zeros(3)))

# Specify that only parameter 3 (γ) should be fitted, keeping S and k fixed (indices 1,2)
sol = solve(bvp, MIRK4(), dt=0.05, abstol=1e-8, reltol=1e-8, maxiters=1000)

# Access the fitted γ value
println("Fitted γ = ", sol.prob.p[3])

sol = solve(bvp, MIRK4(), dt=0.05, abstol=1e-8, reltol=1e-8, maxiters=1000)

# Plot results
plot(sol.t, sol[1,:], label=L"\hat\psi(x)", linewidth=2)
plot!(sol.t, sol[3,:], label=L"\hat\phi(x)", linewidth=2)
xlabel!(L"x")
ylabel!("Solution")
title!("BVP Solution")

sol.prob.p[3]