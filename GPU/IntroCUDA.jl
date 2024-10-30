using CUDA
using BenchmarkTools

const a = Float32(3.14159) # some number (π)
const dim = 1_000_000_000 # 1 billion

# SAXPY on CPU
x = ones(Float32, dim)
y = ones(Float32, dim)
z = zeros(Float32, dim)

typeof(x)

# Perform SAXPY
@btime z .= a .* x .+ y

cpu_time = 499.401 # ms

# SAXPY on GPU
X1 = CUDA.ones(Float32, dim)
Y1 = CUDA.ones(Float32, dim)
Z1 = CUDA.zeros(Float32, dim)

typeof(X1)

# Perform SAXPY
@btime CUDA.@sync Z1 .= a .* X1 .+ Y1

gpu_time = 202.482 # ms

# compare performance
performance = cpu_time / gpu_time

println("My GPU is about ", round(performance, digits=1), " times faster than my CPU for this computation")