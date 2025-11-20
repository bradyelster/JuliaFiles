using BenchmarkTools

#=
function fib(n::Integer)
    if n <= 1
        return 1
    end
    return fib(n - 1) + fib(n - 2)
end
=#

function fibval(::Val{N}) where N
    if N ≤ 1
        1
    else
        fibval(Val{N - 1}()) + fibval(Val{N - 2}())
    end
end

@btime(fibval(Val(46)))