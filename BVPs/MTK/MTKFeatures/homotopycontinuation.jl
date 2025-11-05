using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

# requires dependencies shown here: https://github.com/SciML/ModelingToolkit.jl/blob/8d6be1e304598e12f1e3bf804d051b62cdb11640/test/extensions/homotopy_continuation.jl#L9

# homotopy continuation is an algorithm which turns nonlinear systems
# into polynomials and finds all of their roots

@variables x = 0.25 y = 0.125

eqns = [0 ~ (x^2 - 5x*y + 6y^2) / (x - 0.25)
        0 ~ (x - 0.25) * (x - 0.5)]

@named sys = System(eqns, [x, y], [])
sys = complete(sys)
prob = HomotopyContinuationProblem(sys, [])
alg = HomotopyContinuationJL{true}()

sol = solve(prob, alg)::EnsembleSolution
getproperty.(sol, :u)