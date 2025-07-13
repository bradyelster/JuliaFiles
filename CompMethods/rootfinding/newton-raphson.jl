function newton_raphson(f, f_prime, x0, tol=1e-6, max_iter=100)
    # The Newton-Raphson method requires an initial guess x0 and the function's derivative
    # f_prime. Unlike bisection or false position, it doesn't need an interval.
    
    # Start with initial guess
    x = x0
    
    for iter in 1:max_iter
        # Evaluate function and its derivative at current point
        f_x = f(x)
        fp_x = f_prime(x)
        
        # Check if derivative is too close to zero to avoid division problems
        if abs(fp_x) < 1e-10
            error("Derivative too close to zero. The method may have encountered a local extremum.")
        end
        
        # Calculate next approximation using the Newton-Raphson formula:
        # x_{n+1} = x_n - f(x_n)/f'(x_n)
        # This comes from the linear approximation (tangent line) at the current point
        x_new = x - f_x/fp_x
        
        # Check for convergence
        if abs(x_new - x) < tol
            return x_new
        end
        
        # Update for next iteration
        x = x_new
    end
    
    # If we haven't converged after max_iter iterations, return best approximation
    # but warn the user
    @warn "Maximum iterations reached. The method may not have converged."
    return x
end

# Define our function f(x) = x^2 - 2 and its derivative f'(x) = 2x
f(x) = x^2 - 2
f_prime(x) = 2x

# Find the root starting with initial guess x0 = 2.0
root = newton_raphson(f, f_prime, 2.0)

# Print results
println("Approximate root: ", root)
println("f(root) ≈ ", f(root))

# Let's also see how quickly it converges
function demonstrate_convergence(f, f_prime, x0, tol=1e-6, max_iter=10)
    x = x0
    println("Starting from x0 = ", x0)
    
    for iter in 1:max_iter
        f_x = f(x)
        fp_x = f_prime(x)
        x_new = x - f_x/fp_x
        
        println("Iteration ", iter, ": x = ", x_new, ", f(x) = ", f(x_new))
        
        if abs(x_new - x) < tol
            println("Converged!")
            break
        end
        x = x_new
    end
end

println("\nConvergence demonstration:")
demonstrate_convergence(f, f_prime, 2.0)