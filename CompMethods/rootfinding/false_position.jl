function false_position(f, a, b, tol=1e-6, max_iter=100)
    # Check if the initial interval contains a root
    if f(a) * f(b) >= 0
        error("Function must have opposite signs at interval endpoints")
    end
    
    # Initialize our interval endpoints and their function values
    left, right = a, b
    f_left, f_right = f(left), f(right)
    
    println("\nStarting false position method with interval [", left, ", ", right, "]")
    println("f(left) = ", f_left)
    println("f(right) = ", f_right)
    
    # Keep track of previous approximation for convergence check
    prev_x = left
    
    for iter in 1:max_iter
        # Calculate the next approximation using linear interpolation
        x = left - f_left * (right - left)/(f_right - f_left)
        f_x = f(x)
        
        # Print the current iteration's details
        println("\nIteration ", iter, ":")
        println("x = ", x)
        println("f(x) = ", f_x)
        println("Current interval: [", left, ", ", right, "]")
        
        # Check if we've converged
        if abs(x - prev_x) < tol
            println("\nConverged! Final approximation:")
            println("x ≈ ", x)
            println("f(x) ≈ ", f_x)
            println("Changed by less than tolerance: ", abs(x - prev_x))
            return x
        end
        
        # Update the interval based on where the sign changes
        if f_left * f_x < 0
            # Root is between left and x
            right = x
            f_right = f_x
            println("Root is in left half - updating right endpoint")
        else
            # Root is between x and right
            left = x
            f_left = f_x
            println("Root is in right half - updating left endpoint")
        end
        
        prev_x = x
    end
    
    @warn "Maximum iterations reached. The method may not have converged."
    return prev_x
end

# Define our function f(x) = x^2 - 2
f(x) = x^2 - 2

# Find the root in interval [1, 2]
println("Finding the root of f(x) = x^2 - 2 in the interval [1, 2]")
println("(Expected root is √2 ≈ 1.4142135623730951)")
root = false_position(f, 1.0, 2.0)

# Print final results
println("\nSummary:")
println("Final approximation: ", root)
println("True error: ", abs(root - sqrt(2)))
println("f(root) = ", f(root))