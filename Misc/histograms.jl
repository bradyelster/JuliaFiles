using CairoMakie  # or GLMakie for interactive plots
using Distributions
using DataFrames
using Colors

# Generate similar data to your R example
setseed!(42)
n = 100
group = repeat([1, 2], inner=n ÷ 2)
x = vcat(rand(Normal(1), n ÷ 2), rand(Normal(5), n ÷ 2))
y = vcat(rand(LogNormal(8), n ÷ 2), rand(LogNormal(6), n ÷ 2))
df = DataFrame(x=x, y=y, group=group)

# Create the figure layout
fig = Figure(resolution=(800, 800))

# Create the main scatter plot axis
ax_main = Axis(fig[2, 1],
    xlabel="x", ylabel="y",
    xticks=LinearTicks(5), yticks=LinearTicks(5))

# Create top histogram axis (for x distribution)
ax_top = Axis(fig[1, 1],
    xlimits=ax_main.xlimits, yaxisposition=:right,
    xticks=Visible = false, yticks=Visible = false,
    xgridvisible=false, ygridvisible=false)

# Create right histogram axis (for y distribution)
ax_right = Axis(fig[2, 2],
    ylimits=ax_main.ylimits, xaxisposition=:top,
    xticks=Visible = false, yticks=Visible = false,
    xgridvisible=false, ygridvisible=false)

# Hide decorations for marginal plots
hidedecorations!(ax_top, grid=false)
hidedecorations!(ax_right, grid=false)

# Create scatter plot
colors = [:steelblue, :orange]
for (i, g) in enumerate(unique(df.group))
    group_data = filter(row -> row.group == g, df)
    scatter!(ax_main, group_data.x, group_data.y,
        color=(colors[i], 0.7), markersize=10, label="Group $g")
end

# Create top histogram (x distribution)
for (i, g) in enumerate(unique(df.group))
    group_data = filter(row -> row.group == g, df)
    hist!(ax_top, group_data.x, bins=20, color=(colors[i], 0.5),
        strokecolor=colors[i], strokewidth=1, normalization=:pdf)
end

# Create right histogram (y distribution)
for (i, g) in enumerate(unique(df.group))
    group_data = filter(row -> row.group == g, df)
    hist!(ax_right, group_data.y, bins=20, color=(colors[i], 0.5),
        strokecolor=colors[i], strokewidth=1, direction=:x, normalization=:pdf)
end

# Add legend
Legend(fig[1, 2], ax_main, framevisible=false)

# Adjust layout spacing
colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)

fig