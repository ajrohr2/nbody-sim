# BasicNBodySim.jl Examples

```@meta
CurrentModule = BasicNBodySim
```

## Simulate the Earth-Moon orbit

We first need to create a vector of the bodies we want to simulate. Let's create a planet with mass `5.9722e24` kg centered at (0,0) and a moon with mass `7.348e22` kg a distance of `3.844e8` meters from the center, with a velocity of `1022` m/s in the positive `y` direction. The colors are to tell `CairoMakie` how to render our bodies.

```julia
bodies = [
    Planet(5.9722e24, [0.0, 0.0], [0.0, 0.0], :blue),
    Moon(7.348e22, [3.844e8, 0.0], [0.0, 1022.0], :gray)
]
```

Now we can simulate a year of time (`31536000` seconds) with 1 hour steps (`365*24` steps). This will give us a stable simulation. We don't need a frame for every hour though, so let's use 365 frames (1 frame per day). We want to save our output video to a file called `earth-moon.mp4`.

```julia
simulate(bodies, 31536000.0, 365*24, 365, "earth-moon.mp4")
```

This creates the following video:

```@raw html
<video width="640" height="360" controls>
  <source src="assets/earth-moon.mp4">
  Your browser does not support the video tag.
</video>
```
