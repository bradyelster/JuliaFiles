# Define a type for vectors
struct VectorGA
    components::Vector{Float64}  # Vector components
end

# Overload string representation for nicer printing
Base.show(io::IO, v::VectorGA) = print(io, "Vector(", v.components, ")")

# Wedge product function
function wedge(v1::VectorGA, v2::VectorGA)
    # Ensure both vectors are in the same space
    if length(v1.components) != length(v2.components)
        error("Vectors must have the same dimension for the wedge product.")
    end

    # Compute the bivector components
    n = length(v1.components)
    bivector = []

    # Iterate over all unique combinations of indices (i < j)
    for i in 1:(n - 1)
        for j in (i + 1):n
            push!(bivector, v1.components[i] * v2.components[j] - v1.components[j] * v2.components[i])
        end
    end
    return bivector
end

# Example usage
v1 = VectorGA([1.0, 0.0, 0.0])  # x-axis vector
v2 = VectorGA([0.0, 1.0, 0.0])  # y-axis vector

∧(x, y) = wedge(x, y)
bivector = v1 ∧ v2
println("Wedge product of $v1 and $v2: $bivector")
