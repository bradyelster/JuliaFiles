using QuadGK, Plots, Printf, LaTeXStrings

integrand(x) = 2 * π * (1 / x) * sqrt(1 + (-1 / x^2)^2)
fp(x) = 2 * π * (1 / x) * (sqrt(1 + 1 / x^4) - 1)

# Define range of upper limits to test
upper_limits = 10:10:100

# Arrays to store results for both integrals
integrals_area = Float64[]
errors_area = Float64[]
integrals_fp = Float64[]
errors_fp = Float64[]

# Compute both integrals for each upper limit
println("Computing integrals...")
println("\nArea integral:")
for b in upper_limits
    integral, error = quadgk(integrand, 1, b)
    push!(integrals_area, integral)
    push!(errors_area, error)
    @printf("Upper limit: %4d, Integral: %.6f, Error: %.2e\n", b, integral, error)
end

println("\nfp integral:")
for b in upper_limits
    integral, error = quadgk(fp, 1, b)
    push!(integrals_fp, integral)
    push!(errors_fp, error)
    @printf("Upper limit: %4d, Integral: %.6f, Error: %.2e\n", b, integral, error)
end

# Create the plot with both curves
p = plot(upper_limits, integrals_area,
    xlabel=L"Gabriel's Horn Length $b$",
    ylabel="Integral Value",
    title="Integral Convergence Comparison",
    label="exact value",
    linewidth=2,
    marker=:circle,
    markersize=3,
    legend=:best,
    grid=true,
    minorgrid=false
)

# Add the second integral to the same plot
plot!(p, upper_limits, integrals_fp,
    label="regularization",
    linewidth=2,
    marker=:square,
    markersize=3,
)

display(p)

# Print summary for both integrals
println("\n" * "="^50)
println("Summary - Area Integral:")
println("Last value (b=$(upper_limits[end])): $(integrals_area[end])")
println("\nSummary - fp Integral:")
println("Last value (b=$(upper_limits[end])): $(integrals_fp[end])")
println("="^50)