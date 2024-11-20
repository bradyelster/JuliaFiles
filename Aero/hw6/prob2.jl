using Plots
using LaTeXStrings
using BoundaryValueDiffEq

mutable struct Params
    domain::Tuple{Float64, Float64}
    ϵ::Float64
    α::Int64
    β::Int64
end

function rightBL!(dv, v, p, x)
    u = v[1]
    du = v[2]
    dv[1] = du
    return dv[2] = (du + 2 * x) / inputs.ϵ
end

function bc!(residual, u, p, x)
    residual[1] = u[begin][1] - inputs.α # the solution at the beginning of the time span should be α
    return residual[2] = u[end][1] - inputs.β # the solution at the end of the time span should be β
end

function main()
    global inputs = Params(
        (0.0, 1.0)::Tuple{Float64, Float64},
        0.01::Float64,
        1::Int64,
        2::Int64
    )

    bvp = BVProblem(rightBL!, bc!, [inputs.α, 0], inputs.domain)
    sol = solve(bvp, MIRK4(), dt = 0.001)
    genPlot(sol)
    return Nothing
end

function genPlot(sol)

    theme(:dao::Symbol)
    ϵ = inputs.ϵ
    α = inputs.α
    β = inputs.β

    myplot = plot(
        sol,
        idxs = [1],
        line = (3, :dash),
        label = "Exact",
        title = L"Solution, $\varepsilon = $" * "$ϵ" * L", $\alpha = $" * "$α " * L", $\beta= $" * "$β",
        xlabel = L"$x$",
        ylabel = L"$y(x)$",
        legend = :best,
        titlefontsize = 20,
        tickfontsize = 12,
        legendfontsize = 10,
        yguidefontsize = 15,
        xguidefontsize = 15,
    )
    savefig(myplot, "hw6/Aero_HW6_Prob2.pdf")
    return myplot

end

main()
