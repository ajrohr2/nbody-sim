using BenchmarkTools
import BasicNBodySim

const SUITE = BenchmarkGroup()

SUITE["simulate"] = BenchmarkGroup()

for m in 1:1:5
    bodies = Vector{BasicNBodySim.Body}([BasicNBodySim.Planet(5.9722e24, [0.0, 0.0], [0.0, 0.0], :blue)])
    append!(bodies, [BasicNBodySim.Moon(7.348e22, (2.66e-6).*[-(3.84e8)*sin((2*pi*n)/m), (3.84e8)*cos((2*pi*n/m))], [(3.84e8)*cos((2*pi*n)/m), (3.84e8)*sin((2*pi*n/m))], :black) for n in 1:m])

    SUITE["simulate"][m] = @benchmarkable(BasicNBodySim.simulate(bs, 31536000.0, 365*24*6),evals=1,seconds=60,samples=100,setup=(bs = $bodies))
end
