f(x) = (exp(-x^2)*(x^7 - 5*x^6 + 4*x^4 + x^3- 2*x^2 + 10*x + 6))/(x^9 + 14*x^8 + 19*x^7 - 29*x^6 + 15*x^4 - 13*x^3 + 11*x^2 - 17*x + 5)
N(x) = (x^7 - 5*x^6 + 4*x^4 + x^3- 2*x^2 + 10*x + 6) # polynomial part of the numerator of f(x)

# bisection method algorithm
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

# the function N(x) has roots that are also roots of f(x). The following code shows that
# the bisection method on the same bracketing interval produces the same roots using either function

bisection(N, [-0.7, 0.0], 1e-4)
bisection(N, [1.2, 2.5], 1e-4)
bisection(N, [4, 7], 1e-4)

bisection(f, [-0.7, 0.0], 1e-4)
bisection(f, [1.2, 2.5], 1e-4)
bisection(f, [4, 7], 1e-4)


### Here begins Part 2 ###

# define a function h(x), which is just f(x) for part 2
h2(x) = -x + 2*log(x)+sin(x)
bisection(h2, [1.0, 1.4], 1e-4)
bisection(h2, [2.2, 2.8], 1e-4)