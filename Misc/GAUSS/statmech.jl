using Plots, LaTeXStrings, Printf

xspan = (0, 3)
xvals = LinRange(xspan[1], xspan[2], 1000)
f(x) = x^(-2) * exp(1 / x) / (exp(1 / x) - 1)^2

theme(:dao)
p = plot(
    xvals,
    f.(xvals),
    line=(3, :solid),
    #title=L"$Q^2 = $%$Q2, $Fr = $%$Fr",
    #labels=[L"$\tilde$" L"$$"]
    # title = "Rotating Bead",
    #xlabel=L"$t$",
    #ylabel=L"$\tilde{\phi \ }(\tilde{t})$",
    legend=false,
    #xlabel=L"$\tilde{x}$",
    titlefontsize=20,
    tickfontsize=12,
    legendfontsize=15,
    yguidefontsize=15,
    xguidefontsize=15,
    right_margin=2 * Plots.mm,
    xlims=xspan,
    dpi=300
)
# Add vertical line at t = 1
# p = vline!(p, [1], color=:black, linestyle=:dash, label=false)
#f(t) = 2 * acot(exp(t))
# g(t) = (0.5 * π) * exp(-t)
#tvals = LinRange(0, 20.0, 1000)
# p = plot!(tvals, g.(tvals), line=(2, :dash))
display(p)
savefig("heatcapT.png")