using LaTeXStrings

function riemann(f::Function, endpoints::Tuple, npoints::Int64)
    a, b = endpoints
    Δx = abs(b - a) / npoints
    points = LinRange(a, b, npoints)
    sum = 0
    for i in points
        sum += f(i)*Δx
    end
    return sum
end

f(x) = x^2

function main()
    bounds = (0, 1)
    num_points = 10_000
    sum = riemann(f, bounds, num_points)

    exact = 1//3
    println("The integral from $(bounds[1]) to $(bounds[2]) of x^2 is: $sum")
    println("The error is: ", round(abs(exact - sum)*100, digits=4), "%")
end

main()