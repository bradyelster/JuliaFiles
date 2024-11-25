using Plots
using LaTeXStrings
using Printf
using BoundaryValueDiffEq

struct BVPParams
    domain::Tuple{Float64, Float64}
    ϵ::Vector{Float64}  # Can hold multiple values of epsilon
    α::Int
    β::Int
end

# Differential equation definition
function prob2!(dv, v, x, ϵ)
    @views dv[1] = v[2]
    @views dv[2] = (v[2] + 2x) / ϵ  # Solve for the current epsilon
    return dv
end

# Boundary conditions
function boundary_conditions!(residual, u, params::BVPParams, x)
    residual[1] = u[1][1] - params.α
    residual[2] = u[end][1] - params.β
    return residual
end

# Solve for each epsilon
function solve_bvp(params::BVPParams; dt = 0.01)
    solutions = []  # Array to store solutions
    for ϵ in params.ϵ
        # Define the problem and solve for the current epsilon
        initial_guess = [params.α, 0.0]
        bvp = BVProblem(
            (dv, v, p, x) -> prob2!(dv, v, x, ϵ),  # Closure for the current epsilon
            boundary_conditions!,
            initial_guess,
            params.domain,
            params
        )
        sol = solve(bvp, MIRK4(); dt)
        push!(solutions, (ϵ, sol))  # Store (epsilon, solution) pair
    end
    return solutions
end

# Plot the solutions
function plot_bvp(solutions)
    # Initialize an empty plot
    theme(:dao::Symbol)
    plt = plot(
        xlabel = L"$x$",  # Explicitly set the x-axis label
        ylabel = L"$y(x)$",
        title = "Problem 2 Solution",
        legend = :best,
        titlefontsize = 20,
        tickfontsize = 12,
        legendfontsize = 10,
        yguidefontsize = 15,
        xguidefontsize = 15,
        right_margin = 5 * Plots.mm
    )

    # Loop over solutions and add to plot
    for (ϵ, sol) in solutions
        plot!(
            plt,
            sol,
            idxs = [1],  # Only plot y(x) not y'(x)
            label = L"Exact ($\varepsilon = %$ϵ$)",
            line = (3, :dash),
            xlabel = L"$x$"  # Ensure the x-axis label is overridden for this plot
        )
    end
    return plt
end


# Main function
function main()
    domain = (0.0, 1.0)
    α = 1
    β = 2
    ϵ1 = 0.1
    ϵ2 = 0.01
    params = BVPParams(domain, [ϵ1, ϵ2], α, β)
    solutions = solve_bvp(params)

    points = collect(LinRange(domain[1], domain[2], 1000))
    exact(x, ϵ) = β .+ (β.-α.+1).*(exp.(x ./ ϵ) .- 1) .- x.^2
    asymp1 = exact(points, 0.1)
    asymp2 = exact(points, 0.01)


    fullplot = plot_bvp(solutions)
    plot!(fullplot, points, [asymp1, asymp2])

    # Save the combined plot
    savefig(fullplot, "hw6/Aero_HW6_Prob2.pdf")
end

main()
