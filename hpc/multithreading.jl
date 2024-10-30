using Base.Threads
using BenchmarkTools

nthreads()

a = zeros(12)

Threads.@threads for i = 1:nthreads()
    a[i] = Threads.threadid()
end

a

function sum_single(a)
    # start sum at 0
    s = 0
    for i in a
        # increase the sum by each element's value
        s += i
    end
    # print the total sum
    s
end

function sum_multi(a)
    # partition the threads amongst the elements of a
    chunks = Iterators.partition(a, length(a) ÷ Threads.nthreads())
    # each thread sums its chunk in serial
    tasks = map(chunks) do chunk
        Threads.@spawn sum_single(chunk)
    end
    # place all chunk sums in a vector called chunk_sums
    chunk_sums = fetch.(tasks)
    # return the total sum of each element of chunk_sums
    return sum_single(chunk_sums)
end

@btime sum_single(1:10_000_000)
# 2.260 ns

@btime sum_multi(1:10_000_000)
# 4.308 μs

