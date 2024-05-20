using Profile

function myfunc()
    A = rand(200, 200, 400)
    maximum(A)
end

myfunc()

@profile myfunc()
Profile.print()
@profview myfunc()

using FlameGraphs
data, ldict = Profile.retrieve()
FlameGraphs.save("/benchmarks/experiments/profiling/mittleffsum2_bottleneck.jlprof",Profile.retrieve()...)
