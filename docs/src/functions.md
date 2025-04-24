# BasicNBodySim.jl Simulation Functions

```@meta
CurrentModule = BasicNBodySim
```

## Simulate

```@docs
simulate(all_bodies::Vector{Body}, time::Float64, steps::Int64)
```

## Update
```@docs
update!(current_body::Body, all_bodies::Vector{Body}, dtime::Float64)
```

## Get axis length
```@docs
get_axlen(all_bodies::Vector{Body})
```
