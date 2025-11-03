using LinearAlgebra, SpecialFunctions, Optim, Plots, Statistics

function solve_helmholtz_coil_analog()
    # Parameters
    R = 1.0  # sphere radius
    V = 1.0  # potential magnitude
    L_max = 30  # increased for better accuracy
    
    # Objective function to minimize field variation
    function field_variation(θ0)
        A = calculate_coefficients(θ0[1], R, V, L_max)
        field_std = evaluate_field_uniformity(A, R)
        println("θ0 = $(rad2deg(θ0[1]))°, variation = $field_std")  # Debug output
        return field_std
    end
    
    # Use a different optimization approach - scan the parameter space
    θ0_range = deg2rad.(range(40, 60, length=50))
    variations = Float64[]
    
    for θ0 in θ0_range
        variation = field_variation([θ0])
        push!(variations, variation)
    end
    
    min_idx = argmin(variations)
    θ0_opt = θ0_range[min_idx]
    min_variation = variations[min_idx]
    
    return rad2deg(θ0_opt), min_variation, θ0_range, variations
end

function calculate_coefficients(θ0, R, V, L_max)
    A = zeros(L_max+1)
    cos_θ0 = cos(θ0)
    
    for l in 0:L_max
        if l == 0
            A[1] = 0.0
        else
            # Use orthogonality: A_l = (2l+1)/(2R^l) * ∫V(θ)P_l(cosθ)sinθ dθ
            n_points = 1000
            integral = 0.0
            for i in 1:n_points
                θ = π * (i-1) / (n_points-1)
                x = cos(θ)
                P_l = legendre_poly(x, l)
                
                # Define potential function
                if θ < θ0
                    V_θ = V  # Upper segment
                elseif θ > π - θ0
                    V_θ = -V  # Lower segment  
                else
                    V_θ = 0.0  # Center segment (grounded)
                end
                
                integral += V_θ * P_l * sin(θ)
            end
            integral *= π / (n_points-1)  # Trapezoidal rule normalization
            A[l+1] = (2l+1) / (2 * R^l) * integral
        end
    end
    return A
end

function legendre_poly(x, l)
    if l == 0
        return 1.0
    elseif l == 1
        return x
    else
        P_prev = 1.0
        P_curr = x
        for n in 2:l
            P_next = ((2n-1) * x * P_curr - (n-1) * P_prev) / n
            P_prev, P_curr = P_curr, P_next
        end
        return P_curr
    end
end

function evaluate_field_uniformity(A, R)
    # Sample points on a circle near center
    r_test = 0.01 * R  # Very close to center
    n_points = 36
    field_magnitudes = Float64[]
    
    for i in 1:n_points
        θ = 2π * (i-1) / n_points
        E_r, E_θ = calculate_electric_field(A, r_test, θ, R)
        E_mag = sqrt(E_r^2 + E_θ^2)
        push!(field_magnitudes, E_mag)
    end
    
    field_mean = mean(field_magnitudes)
    field_std = std(field_magnitudes)
    
    # Avoid division by zero
    if field_mean < 1e-12
        return 1.0
    else
        return field_std / field_mean
    end
end

function calculate_electric_field(A, r, θ, R)
    E_r = 0.0
    E_θ = 0.0
    x = cos(θ)
    
    for l in 1:length(A)-1
        if abs(A[l+1]) > 1e-12
            P_l = legendre_poly(x, l)
            
            # Numerical derivative for stability
            Δθ = 1e-6
            P_l_plus = legendre_poly(cos(θ + Δθ), l)
            P_l_minus = legendre_poly(cos(θ - Δθ), l)
            dP_dθ = (P_l_plus - P_l_minus) / (2Δθ)
            
            E_r -= l * A[l+1] * r^(l-1) * P_l
            E_θ -= A[l+1] * r^(l-1) * dP_dθ / r  # Note: 1/r factor included
        end
    end
    
    return E_r, E_θ
end

function visualize_solution(θ0_opt)
    R = 1.0
    V = 1.0
    L_max = 30
    
    # Calculate coefficients for optimal angle
    A_opt = calculate_coefficients(deg2rad(θ0_opt), R, V, L_max)
    
    # 1. Plot potential distribution on sphere
    θ_range = range(0, π, length=100)
    V_θ = [if θ < deg2rad(θ0_opt)
              V
          elseif θ > π - deg2rad(θ0_opt)
              -V
          else
              0.0
          end for θ in θ_range]
    
    p1 = plot(rad2deg.(θ_range), V_θ, 
              xlabel="θ (degrees)", ylabel="Potential V(θ)",
              title="Potential Distribution on Sphere (θ₀ = $(round(θ0_opt, digits=2))°)",
              legend=false, linewidth=2)
    
    # Add vertical lines showing segment boundaries
    vline!([θ0_opt, 180 - θ0_opt], linestyle=:dash, color=:red, label="Segment boundaries")
    
    # 2. Plot field uniformity scan
    θ0_test = range(30, 70, length=50)
    variations = Float64[]
    
    for θ0 in θ0_test
        A = calculate_coefficients(deg2rad(θ0), R, V, L_max)
        variation = evaluate_field_uniformity(A, R)
        push!(variations, variation)
    end
    
    p2 = plot(θ0_test, variations,
              xlabel="θ₀ (degrees)", ylabel="Relative Field Variation",
              title="Field Uniformity vs Segment Angle",
              legend=false, linewidth=2)
    scatter!([θ0_opt], [minimum(variations)], markersize=8, label="Optimum")
    
    # 3. Plot electric field magnitude around center for optimal configuration
    r_plot = 0.1
    ϕ_range = range(0, 2π, length=100)
    E_magnitudes = Float64[]
    
    for ϕ in ϕ_range
        E_r, E_θ = calculate_electric_field(A_opt, r_plot, ϕ, R)
        push!(E_magnitudes, sqrt(E_r^2 + E_θ^2))
    end
    
    p3 = plot(rad2deg.(ϕ_range), E_magnitudes,
              xlabel="Azimuthal Angle ϕ (degrees)", ylabel="|E|",
              title="Electric Field Magnitude at r=$r_plot",
              legend=false, linewidth=2)
    
    # Compare with non-optimal configuration (θ₀ = 45°)
    A_nonopt = calculate_coefficients(deg2rad(45.0), R, V, L_max)
    E_magnitudes_nonopt = Float64[]
    
    for ϕ in ϕ_range
        E_r, E_θ = calculate_electric_field(A_nonopt, r_plot, ϕ, R)
        push!(E_magnitudes_nonopt, sqrt(E_r^2 + E_θ^2))
    end
    
    plot!(p3, rad2deg.(ϕ_range), E_magnitudes_nonopt, 
          linewidth=2, linestyle=:dash, label="θ₀ = 45°")
    
    # 4. Show coefficient spectrum
    l_range = 0:L_max
    p4 = plot(l_range, abs.(A_opt), 
              xlabel="Legendre order l", ylabel="|A_l|",
              title="Coefficient Spectrum",
              yscale=:log10, marker=:circle, legend=false)
    
    # Save individual plots
    savefig(p1, "potential_distribution.png")
    savefig(p2, "field_uniformity_scan.png") 
    savefig(p3, "field_magnitude_comparison.png")
    savefig(p4, "coefficient_spectrum.png")
    
    # Save combined plot
    combined_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))
    savefig(combined_plot, "helmholtz_coil_analysis.png")
    
    println("Plots saved as:")
    println("- potential_distribution.png")
    println("- field_uniformity_scan.png")
    println("- field_magnitude_comparison.png") 
    println("- coefficient_spectrum.png")
    println("- helmholtz_coil_analysis.png (combined)")
    
    return combined_plot
end

function main()
    println("Solving electrostatic Helmholtz coil analog...")
    
    # Run optimization
    println("Running parameter scan...")
    θ0_opt, min_variation, θ0_range, variations = solve_helmholtz_coil_analog()
    
    println("\n=== RESULTS ===")
    println("Optimal angle θ₀ = $(round(θ0_opt, digits=2))°")
    println("Minimum relative field variation: $(round(min_variation, digits=6))")
    
    analytical_guess = rad2deg(acos(1/sqrt(3)))
    println("Analytical 'magic angle': $(round(analytical_guess, digits=2))°")
    println("Difference: $(round(abs(θ0_opt - analytical_guess), digits=2))°")
    
    # Show why we might be getting different results
    println("\n=== ANALYSIS ===")
    println("The optimal angle found is different from the magic angle (54.74°)")
    println("This could be because:")
    println("1. The three-segment configuration behaves differently than the classic two-coil Helmholtz")
    println("2. Numerical precision issues in the coefficient calculation")
    println("3. The 'constant field' criterion might be different for this geometry")
    
    # Visualize the solution and save plots
    println("\nGenerating and saving visualizations...")
    visualize_solution(θ0_opt)
    
    return θ0_opt
end

# Execute
println("Starting computation...")
optimal_theta = main()
println("Computation complete!")

# Additional diagnostic: Check field at center for different angles
println("\n=== FIELD AT CENTER DIAGNOSTIC ===")
R = 1.0
V = 1.0
L_max = 30

test_angles = [40, 45, 50, 54.74, 60]
for θ0 in test_angles
    A = calculate_coefficients(deg2rad(θ0), R, V, L_max)
    variation = evaluate_field_uniformity(A, R)
    # Calculate field at center (should be zero for symmetric configurations)
    E_r_center, E_θ_center = calculate_electric_field(A, 0.001, 0.0, R)
    println("θ₀ = $(θ0)°: variation = $(round(variation, digits=4)), E_center = $(round(E_r_center, digits=8))")
end

# Create a detailed analysis plot showing the variation curve
println("\n=== CREATING DETAILED ANALYSIS PLOT ===")
detailed_θ0_range = range(30, 70, length=100)
detailed_variations = Float64[]

for θ0 in detailed_θ0_range
    A = calculate_coefficients(deg2rad(θ0), 1.0, 1.0, 30)
    variation = evaluate_field_uniformity(A, 1.0)
    push!(detailed_variations, variation)
end

detailed_plot = plot(detailed_θ0_range, detailed_variations,
                     xlabel="θ₀ (degrees)", ylabel="Relative Field Variation",
                     title="Detailed Field Uniformity Analysis",
                     linewidth=2, legend=false)

# Mark key points
min_idx = argmin(detailed_variations)
θ0_min = detailed_θ0_range[min_idx]
scatter!([θ0_min], [detailed_variations[min_idx]], 
         markersize=8, label="Numerical optimum: $(round(θ0_min, digits=2))°")

# Mark magic angle
magic_angle = rad2deg(acos(1/sqrt(3)))
scatter!([magic_angle], [detailed_variations[findfirst(x -> x >= magic_angle, detailed_θ0_range)]],
         markersize=8, label="Magic angle: $(round(magic_angle, digits=2))°")

savefig(detailed_plot, "detailed_analysis.png")
println("Detailed analysis plot saved as: detailed_analysis.png")

println("\nAll plots have been saved to PNG files in the current directory.")