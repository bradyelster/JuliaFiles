using BoundaryValueDiffEq
using Plots
using LaTeXStrings

L = 10;
tspan = (0, L);

function schrodinger!(dψ, ψ, p, t)
    dψ[1] = ψ[2]
    dψ[2] = (t^2 - 2 * p[1]) * ψ[1]
end

function bca!(res, ψ, p)
    res[1] = ψ[2]
end

function bcb!(res, ψ, p)
    res[1] = ψ[1]
end

guess(p, t) = [exp(-t^2 / 2); -t * exp(-t^2 / 2)]  # Better initial guess

bvp = TwoPointBVProblem(schrodinger!, (bca!, bcb!), guess, tspan, [6.2],
    bcresid_prototype=(zeros(1), zeros(1)), fit_parameters=true)

sol = solve(bvp, MIRK4(), dt=0.05)

sol.prob.p

plot(sol.t, sol[1, :], label=L"$\psi(x)$", xlabel=L"$x$", ylabel=L"$\psi$",
    xticks=-L:5:L, lw=2)