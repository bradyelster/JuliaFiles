using BoundaryValueDiffEq, Plots, LaTeXStrings

L = 10.0
tspan = (0.0, L)
F(t) = tanh(t)
Fpp(t) = -2 * tanh(t) * sech(t)^2

k = 0.5
S = 100
γ_guess = S^(-2/5)

function tearing!(du, u, p, t)
    ψ, φ, ψp, φp = u
    du[1] = ψp
    du[2] = φp
    du[3] = (S*p[1]+k^2)*ψ - S*p[1]*F(t)*φ
    du[4] = (-F(t)/p[1]^2)*(du[3]-k^2*ψ - Fpp(t)*ψ/F(t)) + k^2*φ
end

function bca!(res, u, p)
    res[1] = u[3]       # ψ'(0)=0
    res[2] = u[2]       # φ(0)=0
end

function bcb!(res, u, p)
    res[1] = u[1]       # ψ(L)=0
    res[2] = u[2]       # φ(L)=0
end

function initial_guess(p, t)
    [exp(-k*t) + F(t)/k; 
    1/t;
    -k * exp(-k*t);
    -1/t^2]
end

bvp = TwoPointBVProblem(
    tearing!,
    (bca!, bcb!),
    initial_guess,
    tspan,
    [γ_guess],
    bcresid_prototype=(zeros(2), zeros(2)),
    fit_parameters=true
)

# Solve
sol = solve(bvp, MIRK3(), dt=0.05)

# print the estimated value of γ which satisfies the BCs
γ = sol.prob.p[1]

plot(sol, idxs=(0, 1), label=L"ψ(t)", continuity=:right)
plot!(sol, idxs=(0, 2), label=L"φ(t)", xlabel="t", legend=:topright, continuity=:right)
