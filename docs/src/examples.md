# BasicNBodySim.jl Examples

```@meta
CurrentModule = BasicNBodySim
```

## Simulate the Earth-Moon orbit

We first need to create a vector of the bodies we want to simulate. Let's create a planet with mass `5.9722e24` kg centered at (0,0) and a moon with mass `7.348e22` kg a distance of `3.844e8` meters from the center, with a velocity of `1022` m/s in the positive `y` direction.

```julia
bodies = [
    Planet(5.9722e24, [0.0, 0.0], [0.0, 0.0], :blue),
    Moon(7.348e22, [3.844e8, 0.0], [0.0, 1022.0], :gray)
]
```

Now we can simulate a year of time (`31536000` seconds) with 1 hour steps (`365*24` steps). This will give us a stable simulation.

```julia
simulate(bodies, 31536000.0, 365*24, 365, "earth-moon.mp4")
```

This creates the following video:
