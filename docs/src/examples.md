# BasicNBodySim.jl Examples

```@meta
CurrentModule = BasicNBodySim
```

## Simulate the Earth-Moon orbit

We first need to create a vector of the bodies we want to simulate. Let's create a planet with mass `5.9722e24` kg centered at (0,0) and a moon with mass `7.348e22` kg a distance of `3.844e8` meters from the center, with a velocity of `1022` m/s in the positive `y` direction.

```julia
bodies = [
    Planet(5.9722e24, [0.0, 0.0], [0.0, 0.0], :blue),
    Moon(7.348e22, [0.0, 1022.0], [3.844e8, 0.0], :gray)
]
```

Now we can simulate a year of time (`31536000` seconds) with 1 hour steps (`365*24` steps). This will give us a stable simulation.

```julia
simulate(bodies, 31536000.0, 365*24)
```

Now if we check the position and velocities of the Planet and Moon objects, we get the following:

```julia
julia> bodies[1].position
2-element Vector{Float64}:
 1.110735831933807e7
 3.9101189854465204e8
```

```julia
julia> bodies[1].velocity
2-element Vector{Float64}
  1.9241461244645808
 24.82682974654261
```

Our Moon certainly moved, however it's difficult to see if the motion looks correct. To do this, we can create an animation using `CairoMakie.jl`, as explained in the next example.

## Save an animation of a simulation

Note: This requires the use of `CairoMakie.jl` to create the animations.

We can now simulate the n-body problem, but it would be nice to have a visual sanity check to determine what the motion looks like. To do this, we will use `CairoMakie.jl`'s `record` function along with our `update!` function. Consider the following function:

```julia
using CairoMakie
using BasicNBodySim

function simulate_with_save(all_bodies::Vector{Body}, time::Float64, steps::Int64, frames::Int64, name::String)
    dtime = time / steps

    positions = [Point2f(p.position) for p in all_bodies]

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=1)

    axlen = get_axlen(all_bodies)

    xlims!(ax, -axlen, axlen)
    ylims!(ax, -axlen, axlen)
    sc = scatter!(ax, positions, color=[b.color for b in all_bodies])

    framerate = 30
    frame_interval = div(steps, frames)

    record(fig, name, 1:frames; framerate = framerate) do frame
        for _ in 1:frame_interval
            for body in all_bodies
                update!(body, all_bodies, dtime)
            end
        end
        central_pos = all_bodies[1].position
        new_positions = [Point2f(p.position .- central_pos) for p in all_bodies]
        sc[1][] = new_positions

        new_axlen = get_axlen(all_bodies)
        if new_axlen > axlen
            axlen = new_axlen
            xlims!(ax, -axlen, axlen)
            ylims!(ax, -axlen, axlen)
        end
    end
end
```

This function uses the `update!` function to step through individual updates of our system. We determine how many frames to record, then update the system and record the positions of each body every `frames` step. Then we save the animation to a file specified in `name`.

Note: it's not necessary to record every frame, as this doesn't present new information visually and creates a very long video. 

If we record the same example we gave before, we get the following video:

```@raw html
<video width="640" height="360" controls>
  <source src="assets/earth-moon.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
```
