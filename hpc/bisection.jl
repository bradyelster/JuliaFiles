function bisection(f::Function, interval::Vector, tol::Float64)
    a, b = interval
    iter = 0
    
    while (b - a) / 2 > tol  # Stop when the interval size is within tolerance
        c = (a + b) / 2      # Midpoint
        if f(a) * f(c) < 0   # If signs are different:
            b = c            # root is in the [a, c] interval
        else                 # else:
            a = c            # root is in the [c, b] interval
        end
        iter += 1
    end
    
    c = (a + b) / 2  # Final midpoint approximation
    return c, iter
end

f(x) = sin(cos(exp(x)))

function main()
    tolerance = 1e-8
    interval = [0, 1]
    root, iterations = bisection(f, interval, tolerance)
    println("root is approximately: $root")
    println("converged in $iterations iterations") 
end

main()