# Hello friend
# it's been a while

using BoundaryValueDiffEq
using Plots
using LaTeXStrings

L = 10;
tspan = (-L, L);

function schrodinger!(dψ, ψ, p, t)
    dψ[1] = ψ[2]
    dψ[2] = (t^2 - 2*p[1]) * ψ[1]
end

# for eigenvalue problems, the solution is only defined up to a multiplicative constant by the boundary conditions
# to resolve this, specifiy another condition like ψ(0) = 1 to determine eigenvalue (E) uniquely
function bca!(res, ψ, p)
    res[1] = ψ[1] 
    # res[2] = ψ[1] - 1.0 # extra constraint that ψ(0)=1
end

function bcb!(res, ψ, p)
    res[1] = ψ[1]
end

guess(p, t) = [cos(π*t / (2*L)); -sin(π*t / (2*L))] # initial guess of [ψ, ψ'], notice this satisfies the BCs exactly

bvp = TwoPointBVProblem(schrodinger!, (bca!, bcb!), guess, tspan, [3.5],
    bcresid_prototype = (zeros(1), zeros(1)), fit_parameters = true)

sol = solve(bvp, MIRK4(), dt = 0.05)

plot(sol.t, sol[1,:], label=L"$\psi(x)$", xticks=-L:5:L, lw=2)

sol.prob.p