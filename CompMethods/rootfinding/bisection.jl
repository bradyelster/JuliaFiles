function bisection(f, a, b, tol=1e-6, max_iter=100)
    # Check if the initial interval has a root
    if f(a) * f(b) >= 0
        error("Function must have opposite signs at interval endpoints")
    end
    
    # Initialize iteration counter and interval endpoints
    left = a
    right = b
    f_left = f(left)
    f_right = f(right)
    
    println("\nStarting bisection method with interval [", left, ", ", right, "]")
    println("f(left) = ", f_left)
    println("f(right) = ", f_right)
    
    # Track the width of our interval for convergence checking
    interval_width = right - left
    
    for iter in 1:max_iter
        # Calculate interval information
        interval_width = right - left
        mid = (left + right) / 2
        f_mid = f(mid)
        
        # Print current iteration details
        println("\nIteration ", iter, ":")
        println("Current interval: [", left, ", ", right, "]")
        println("Interval width: ", interval_width)
        println("Midpoint = ", mid)
        println("f(midpoint) = ", f_mid)
        
        # Check if we found exact root (unlikely with floating-point)
        if f_mid == 0
            println("\nExact root found!")
            return mid
        end
        
        # Determine which half of interval contains the root
        if f_left * f_mid < 0
            # Root is in left half
            right = mid
            println("Root is in left half - updating right endpoint")
        else
            # Root is in right half
            left = mid
            f_left = f_mid
            println("Root is in right half - updating left endpoint")
        end
        
        # Check for convergence based on interval width
        if interval_width < tol
            println("\nConverged! Final interval width ", interval_width, " is less than tolerance ", tol)
            final_approximation = (left + right) / 2
            println("Final approximation: ", final_approximation)
            println("f(approximation) = ", f(final_approximation))
            return final_approximation
        end
    end
    
    # If we reach here, we've hit max iterations without converging
    @warn "Maximum iterations reached. The method may not have converged."
    return (left + right) / 2
end

# Define our function f(x) = x^2 - 2
f(x) = x^2 - 2

# Find the root in interval [1, 2]
println("Finding the root of f(x) = x^2 - 2 in the interval [1, 2]")
println("(Expected root is √2 ≈ 1.4142135623730951)")
root = bisection(f, 1.0, 2.0)

# Print final results
println("\nSummary:")
println("Final approximation: ", root)
println("True error: ", abs(root - sqrt(2)))
println("f(root) = ", f(root))