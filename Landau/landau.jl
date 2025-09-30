using DifferentialEquations
using Makie
using CairoMakie  # Using CairoMakie as requested
using Colors, ColorSchemes
using StatsBase
using Random
using SpecialFunctions: erf
using LinearAlgebra
using LaTeXStrings  # For better fraction formatting

# Constants
const k = 2π  # three full eyes in the simulation domain
const ω = 4π  # two full oscillations in 1s
const g0 = 0.8  # wave amplitude
const t_on = 1.0  # amount of time before wave appears
const Δt_on = 0.07  # amount of time it takes wave to ramp up
const v_thermal = 0.8 * ω / k  # ensure a high gradient at the wave velocity

# Grid definitions
const x_grid = range(-1.5, 1.5, length=361)  # normalized spatial coordinates
const v_grid = range(-0.6 * ω / k, 2.1 * ω / k, length=201)  # velocity bounds

const num_samples = 1_000_000
const frame_rate = 24  # Increased for smoother animation
const duration = t_on + 2.0

# Custom colormap similar to the Python version
const colormap = cgrad(:thermal, rev=true)

function periodicize(x, minimum, maximum)
    return minimum .+ mod.(x .- minimum, maximum - minimum)
end

function format_as_fraction(coefficient, numerator, denominator)
    isapprox(coefficient, 0, atol=1e-6) && return L"0"
    coefficient < 0 && return L"-%$(format_as_fraction(abs(coefficient), numerator, denominator))"

    # Check for common fractions
    common_fractions = [
        (0.5, L"\frac{1}{2}"),
        (0.25, L"\frac{1}{4}"),
        (0.75, L"\frac{3}{4}"),
        (0.333, L"\frac{1}{3}"),
        (0.666, L"\frac{2}{3}"),
        (0.2, L"\frac{1}{5}"),
        (0.4, L"\frac{2}{5}"),
        (0.6, L"\frac{3}{5}"),
        (0.8, L"\frac{4}{5}")
    ]

    for (val, str) in common_fractions
        if isapprox(coefficient, val, atol=1e-3)
            return L"%$str \frac{%$numerator}{%$denominator}"
        end
    end

    # Fallback to decimal representation
    return L"%$(round(coefficient, digits=2)) \frac{%$numerator}{%$denominator}"
end

function derivative!(du, u, p, t)
    x = @view u[1:2:end]
    v = @view u[2:2:end]
    dxdt = @view du[1:2:end]
    dvdt = @view du[2:2:end]

    field_on = p[1]

    @. dxdt = v
    if field_on
        @. dvdt = g0 * (1 + erf((t - t_on) / Δt_on)) / 2 * sin(k * x - ω * t)
    else
        fill!(dvdt, 0.0)
    end
end

function generate_particles()
    Random.seed!(1)
    v0 = randn(num_samples) * v_thermal  # maxwellian initial distribution

    # Exclude particles off screen
    in_bounds = (v0 .> v_grid[1] - 0.1 * ω / k) .& (v0 .< v_grid[end] + 0.1 * ω / k)
    v0 = v0[in_bounds]
    x0 = rand(length(v0)) * (x_grid[end] - x_grid[1]) .+ x_grid[1]

    return x0, v0
end

function solve_equations(x0, v0, field_on=true)
    u0 = zeros(2 * length(x0))
    u0[1:2:end] .= x0
    u0[2:2:end] .= v0

    tspan = (0.0, duration)
    t = range(tspan..., step=1 / frame_rate)

    prob = ODEProblem(derivative!, u0, tspan, [field_on])
    sol = solve(prob, Tsit5(), saveat=t)

    # Extract positions and velocities
    x = sol[1:2:end, :]
    v = sol[2:2:end, :]

    return x, v, sol.t
end

function create_phase_space_plot(x_grid, v_grid, t, x, v, field_on, wave_frame, trajectories, i)
    fig = Figure(size=(1200, 900))  # Fixed deprecation warning

    # Create layout
    grid = fig[1, 1] = GridLayout()

    # Top row: Potential and Field
    ax_V = Axis(grid[1, 1], ylabel=L"Potential", ylabelsize=20, ylabelcolor=:orange)
    ax_E = Axis(grid[1, 1], ylabel=L"Field", ylabelsize=20, ylabelcolor=:purple,
        yaxisposition=:right)

    # Right column: Velocity distribution
    ax_v = Axis(grid[2, 2], xlabel=L"Distribution", xlabelsize=20, xlabelcolor=RGB(0.38, 0.56, 0.62))

    # Main plot: Phase space
    ax_image = Axis(grid[2, 1], xlabel=L"Position", ylabel=L"Velocity", ylabelsize=20)

    # Time display
    textbox = Label(grid[1, 2], L"t = %$(round(t[i], digits=1))\,s",
        halign=:left, valign=:top, fontsize=20)

    # Move with the wave (or not)
    x_plot = wave_frame ? x_grid .+ ω / k * t[i] : x_grid

    # Calculate field and potential
    E = field_on ? g0 * (1 + erf((t[i] - t_on) / Δt_on)) / 2 * sin.(k * x_plot .- ω * t[i]) : zeros(length(x_plot))
    ϕ = field_on ? g0 * (1 + erf((t[i] - t_on) / Δt_on)) / 2 * cos.(k * x_plot .- ω * t[i]) / k : zeros(length(x_plot))

    # Plot field and potential
    lines!(ax_E, x_plot, E, color=:purple, linewidth=1.4)
    lines!(ax_V, x_plot, ϕ, color=:orange, linewidth=1.4, linestyle=:dot)

    # Set limits
    ylims!(ax_E, -1.1 * g0, 1.1 * g0)
    ylims!(ax_V, -1.4 * g0 / k, 1.4 * g0 / k)

    # Plot velocity distribution
    f_v = fit(Histogram, v[:, i], v_grid[1:2:end])
    barplot!(ax_v, f_v.weights ./ diff(f_v.edges[1]), f_v.edges[1][1:end-1],
        direction=:x, color=RGB(0.38, 0.56, 0.62))

    # Theoretical Maxwellian
    theoretical = length(v[:, i]) / sqrt(2π) / v_thermal * exp.(-(v_grid ./ v_thermal) .^ 2 / 2)
    lines!(ax_v, theoretical, v_grid, color=:black, linewidth=1.0, linestyle=:dash)

    xlims!(ax_v, 0, 0.43 * num_samples / v_thermal)

    # Plot 2D distribution function
    r_particle = (x_grid[2] - x_grid[1]) * 1.5
    conv_points = 7
    dx_values = range(-r_particle, r_particle, length=conv_points)
    dy_values = range(-r_particle, r_particle, length=conv_points)

    image = zeros(length(x_grid) - 1, length(v_grid) - 1)
    for dx in dx_values, dy in dy_values
        if hypot(dx, dy) < r_particle * 1.1
            dv = dy / (x_grid[2] - x_grid[1]) * (v_grid[2] - v_grid[1])
            x_periodic = periodicize(x[:, i] .+ dx, x_grid[1], x_grid[end])
            h = fit(Histogram, (x_periodic, v[:, i] .+ dv), (x_grid, v_grid))
            image .+= h.weights
        end
    end

    # Find max value for consistent color scaling
    max_val = maximum(image) * 0.75
    max_val *= trajectories ? 1.5 : 1.0

    heatmap!(ax_image, x_grid[1:end-1], v_grid[1:end-1], image',
        colormap=colormap, colorrange=(0, max_val))

    xlims!(ax_image, x_grid[1] + 0.21, x_grid[end] - 0.21)

    # Format axes
    ax_image.xticks = -1.5:0.5:1.5
    if field_on
        ticks = range(v_grid[1], v_grid[end], step=0.5 * ω / k)
        ax_image.yticks = (ticks, [format_as_fraction(t / (ω / k), L"\omega", L"k") for t in ticks])
    else
        ax_image.yticks = range(v_grid[1], v_grid[end], step=1.0)
    end

    # Plot trajectories if enabled
    if trajectories
        X_grid = [x for x in x_grid, v in v_grid]
        V_grid = [v for x in x_grid, v in v_grid]

        V_max_grid = sqrt.((V_grid .- ω / k) .^ 2 .+ 2 * (ϕ .- minimum(ϕ))')
        v_max_separatrix = 2 * sqrt(g0 / k * (1 + erf((t[i] - t_on) / Δt_on)) / 2)
        v_max_contours = vcat(-4/5:2/5:8) .* 2 * sqrt(g0 / k) .+ v_max_separatrix

        if t[i] >= t_on + Δt_on
            contour!(ax_image, x_grid, v_grid, V_max_grid,
                levels=v_max_contours, linewidth=[1.4, 1.4, 1.4, 0.7, 0.7, 0.7],
                color=:black)
        end

        # Add horizontal line for separatrix if not shown
        if !any(V_max_grid .< v_max_contours[3])
            hlines!(ax_image, [ω / k], color=:black, linewidth=0.7)
        end
    elseif wave_frame
        hlines!(ax_image, [ω / k], color=:black, linewidth=1.0, linestyle=:dash)
    end

    # Link axes
    linkxaxes!(ax_image, ax_V, ax_E)
    linkyaxes!(ax_image, ax_v)

    return fig
end

function generate_animation()
    println("Generating particles...")
    x0, v0 = generate_particles()

    println("Solving equations...")
    x, v, t = solve_equations(x0, v0)

    println("Creating animation frames...")
    frames = []

    # First pass to determine max histogram value for consistent scaling
    max_hist_value = 0.0
    sample_indices = round.(Int, range(1, length(t), length=min(10, length(t))))

    for i in sample_indices
        r_particle = (x_grid[2] - x_grid[1]) * 1.5
        conv_points = 7
        dx_values = range(-r_particle, r_particle, length=conv_points)
        dy_values = range(-r_particle, r_particle, length=conv_points)

        image = zeros(length(x_grid) - 1, length(v_grid) - 1)
        for dx in dx_values, dy in dy_values
            if hypot(dx, dy) < r_particle * 1.1
                dv = dy / (x_grid[2] - x_grid[1]) * (v_grid[2] - v_grid[1])
                x_periodic = periodicize(x[:, i] .+ dx, x_grid[1], x_grid[end])
                h = fit(Histogram, (x_periodic, v[:, i] .+ dv), (x_grid, v_grid))
                image .+= h.weights
            end
        end
        max_hist_value = max(max_hist_value, maximum(image))
    end

    vmax_value = max_hist_value * 0.75 * 1.5  # Using 1.5 for trajectories

    # Second pass to create frames
    for i in 1:length(t)
        if i % 10 == 0
            println("  Processing frame $i/$(length(t))...")
        end

        fig = create_phase_space_plot(x_grid, v_grid, t, x, v, true, true, true, i)

        # Save frame to memory
        push!(frames, fig)
    end

    println("Saving animation...")
    save("landau_damping.gif", frames, fps=frame_rate)
    println("Done! Saved 'landau_damping.gif'")
end

# Run the simulation
generate_animation()