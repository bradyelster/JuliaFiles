# ζ(s) ≔ ∑ 1/n^s from n=1 to n = ∞. 
using Plots, LaTeXStrings

function zeta(order::Int, N::Int)
    approx = 0
    for n in 1:N
        approx += 1 / n^order
    end
    return approx
end

# Set the order for zeta function
order = 3

# Generate a range of N values
N_values = 1:500

# Calculate partial sums for each N
partial_sums = [zeta(order, N) for N in N_values]

# The exact value of ζ(4) = π^4/90, ζ(2) = π^2/6
# Uncomment if you want to add exact values for specific orders
# exact_values = Dict(2 => π^2/6, 4 => π^4/90, 6 => π^6/945)
# zeta_exact = get(exact_values, order, nothing)

# Display the final approximation
println("Partial sum at N=$(N_values[end]): ", zeta(order, N_values[end]))
# if !isnothing(zeta_exact)
#     println("Exact value ζ($order): ", zeta_exact)
#     println("Error: ", abs(zeta(order, N_values[end]) - zeta_exact))
# end

# Create the plot
theme(:rose_pine)
plot(N_values, partial_sums,
    label="Partial Sum",
    xlabel=L"N",
    ylabel=L"\sum_{n=1}^{N} \frac{1}{n^{%$order}}",
    title=L"Convergence to $\zeta(%$order)$",
    linewidth=2,
    legend=:bottomright,
    xscale=:log10,
    left_margin=2Plots.mm
)

# Optionally add exact value line
# if !isnothing(zeta_exact)
#     hline!([zeta_exact], 
#            label=L"\zeta(%$order)",
#            linestyle=:dash,
#            linewidth=2,
#            color=:red)
# end