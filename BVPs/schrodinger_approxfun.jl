using ApproxFun

x = Fun(-8 .. 8)
V = x^2 / 2
L = -𝒟^2 / 2 + V
S = space(x)
B = Dirichlet(S)
λ, v = ApproxFun.eigs(B, L, 100, tolerance=1E-5);

import Plots
using LinearAlgebra: norm
p = Plots.plot(V, legend=false, ylim=(-Inf, λ[7]))
for k = 1:8
    Plots.plot!(real(v[k]) + 2.2 * real(λ[k]), xticks=-7:2:7)
end
p