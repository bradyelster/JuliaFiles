using Plots
using LaTeXStrings
using BoundaryValueDiffEq

# Immutable struct for better performance
struct BVPParams
    domain::Tuple{Float64, Float64}
    ϵ::Float64
    α::Int
    β::Int
end

# Use multiple dispatch instead of global variable
function right_bl!(dv, v, params::BVPParams, x)
    @views dv[1] = v[2]
    @views dv[2] = (v[2] + 2x) / params.ϵ
    return dv
end

function boundary_conditions!(residual, u, params::BVPParams, x)
    residual[1] = u[1][1] - params.α
    residual[2] = u[end][1] - params.β
    return residual
end

function generate_plot(sol, params::BVPParams)
    return plot(
        sol,
        idxs = [1],
        line = (3, :dash),
        label = "Exact",
        title = L"Solution, $\varepsilon = %$(params.ϵ)$, $\alpha = %$(params.α)$, $\beta = %$(params.β)$",
        xlabel = L"$x$",
        ylabel = L"$y(x)$",
        legend = :best,
        titlefontsize = 20,
        tickfontsize = 12,
        legendfontsize = 10,
        yguidefontsize = 15,
        xguidefontsize = 15,
        theme = :dao
    )
end

function solve_bvp(params::BVPParams; dt = 0.001)
    initial_guess = [params.α, 0.0]
    bvp = BVProblem(right_bl!, boundary_conditions!, initial_guess, params.domain, params)
    sol = solve(bvp, MIRK4(); dt)

    plot = generate_plot(sol, params)
    savefig(plot, "hw6/Aero_HW6_Prob2_improved.pdf")

    return (solution = sol, plot = plot)
end

# Example usage
params = BVPParams((0.0, 1.0), 0.01, 1, 2)
result = solve_bvp(params)
