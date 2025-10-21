# Rotating Bead on Hoop - plots for GAUSS Talk

using OrdinaryDiffEq, Plots, LaTeXStrings

# Parameters
Fr = 0.1
Q2 = 25

# Define the ODE system
function phi_ode!(du, u, p, t)
    # u[1] = φ, u[2] = φ'
    φ, φdot = u
    du[1] = φdot
    du[2] = -(1 / Q2) * (φdot + (1 - Fr * cos(φ)) * sin(φ))
end

# Initial conditions
u0 = [π / 2, 0.0]  # φ(0) = 1, φ'(0) = 0
tspan = (0.0, 150.0)

# Define and solve the problem
prob = ODEProblem(phi_ode!, u0, tspan)
sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-8)

theme(:dao)
p = plot(
    sol,
    idxs=[1],
    line=(3, :solid),
    title=L"$Q^2 = $%$Q2, $Fr = $%$Fr",
    #labels=[L"$\tilde$" L"$$"]
    # title = "Rotating Bead",
    #xlabel=L"$t$",
    #ylabel=L"$\tilde{\phi \ }(\tilde{t})$",
    legend=false,
    xlabel="",
    titlefontsize=20,
    tickfontsize=12,
    legendfontsize=15,
    yguidefontsize=15,
    xguidefontsize=15,
    right_margin=2 * Plots.mm,
    xlims=tspan,
    #yticks=(
    #[0, (π - 2) / 2, π / 2, π, 2π],
    #[L"0", L"\pi/2 - 1", L"\pi/2", L"\pi", L"2\pi"]),
    dpi=300
)
# Add vertical line at t = 1
# p = vline!(p, [1], color=:black, linestyle=:dash, label=false)
#f(t) = 2 * acot(exp(t))
# g(t) = (0.5 * π) * exp(-t)
#tvals = LinRange(0, 20.0, 1000)
# p = plot!(tvals, g.(tvals), line=(2, :dash))
display(p)
savefig("highQualityFactor.png")


# 4 regimes
# 1. VG: viscous-gravitational (Q << 1, Fr << 1)
# 2. VC: viscous-centrifugal (Q << 1, Fr >> 1)
# 3. IG: inertial-gravitational (Q >> 1, Fr << 1)
# 4. IC: inertial-centrifugal (Q >> 1, Fr >> 1)