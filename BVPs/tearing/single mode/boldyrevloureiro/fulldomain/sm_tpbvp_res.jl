using BoundaryValueDiffEq, Plots, LaTeXStrings

const L = 15.0
tspan = (-L, L)
@inline f(t) = tanh(t)

const k = 0.5
const S = 100
γ_guess = 0.01326

ψouter(t) = exp(-k * abs(t)) * (1 + f(abs(t)) / k)
ϕouter(t) = ψouter(t) / f(t)

function tearing!(du, u, p, t)
    ψ, ϕ, ψ1, ϕ1 = u
    γ = p[1]

    du[1] = ψ1
    du[2] = ϕ1

    # ψ'' equation
    du[3] = (-S * k * f(t)) * ϕ + (S * γ - k^2) * ψ

    # ϕ'' equation
    du[4] = ((k^2 * S * f(t)^2) / γ) * ϕ -
            (k * f(t) / γ) * (-2 * k^2 + S * γ + 2 * (sech(t)^2)) * ψ
end

# boundary conditions: use the interpolant u(t_point) to evaluate fields
function bc!(res, u, p, t)
    res[1] = u(-L)[1] - ψouter(-L)
    res[2] = u(-L)[2] - ϕouter(-L)
    res[3] = u(L)[1] - ψouter(L)
    res[4] = u(L)[2] - ϕouter(L)
    res[5] = u(0.0)[1] - 1.0   # normalization: ψ(x=0) = 1
end

function initial_guess(p, t)
    exp_t2 = exp(-t^2) # guess a Gaussian distribution (decays to zero and ψ(0) = 1)
    t_exp = t * exp_t2
    [exp_t2, t_exp, -2 * t_exp, (1 - 2 * t^2) * exp_t2]
end

bvp = BVProblem(
    BVPFunction(tearing!, bc!; bcresid_prototype=zeros(5)),
    initial_guess, tspan, [γ_guess], fit_parameters=true)

@time sol = solve(bvp, MIRK4(), dt=0.01,
    saveat=0.1,
    adaptive=true,
    progress=true,
    verbose=true
)

println("γ = ", sol.prob.p[1])

plot(sol, idxs=(0, 1), label=L"ψ(t)", lw=2)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:best, lw=2)