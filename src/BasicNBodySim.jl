module BasicNBodySim

using LinearAlgebra

abstract type Body end

mutable struct Planet <: Body
    mass::Float64
    velocity::Vector{Float64}
    position::Vector{Float64}
    color::Symbol

    function Planet(mass::Float64, velocity::Vector{Float64}, position::Vector{Float64}, color)
        new(min(mass, 1.899e27), velocity, position, color)
    end
end

# Some common values:
# Earth: 5.9722e24 kg
# Mars: 6.417e23 kg
# Jupiter: 1.898e27 kg

mutable struct BlackHole <: Body
    mass:: Float64
    velocity::Vector{Float64}
    position::Vector{Float64}
    color::Symbol
    function BlackHole(mass::Float64, velocity::Vector{Float64}, position::Vector{Float64}, color)
        new(min(mass, 3e32), velocity, position, color)
    end
end

mutable struct Moon <: Body
    mass::Float64
    velocity::Vector{Float64}
    position::Vector{Float64}
    color::Symbol
    function Moon(mass::Float64, velocity::Vector{Float64}, position::Vector{Float64}, color)
        new(min(mass, 1.35e23), velocity, position, color)
    end
end

# Some common values:
# Earth's moon: 7.348e22 kg, 1022 m/s
# Titan: 1.35e23 kg, 1.2 million km from Saturn, 5570 m/s

function update!(current_body::Body, all_bodies::Vector{Body}, dtime::Float64)
    surrounding_bodies = filter(e->e != current_body, all_bodies)

    G = 6.6743e-11
    force = [0.0, 0.0]
    for neighbor in surrounding_bodies
        displacement = current_body.position .- neighbor.position
        force .+= -G*(current_body.mass * neighbor.mass)/norm(displacement)^3 * displacement
    end
    acceleration = force ./ current_body.mass
    current_body.position .+= current_body.velocity .* dtime .+ 0.5 .* acceleration .* dtime^2

    new_force = zeros(2)
    for neighbor in surrounding_bodies
        displacement = current_body.position .- neighbor.position
        new_force .+= -G*(current_body.mass * neighbor.mass)/norm(displacement)^3 * displacement
    end
    new_accel = new_force ./ current_body.mass

    current_body.velocity .+= 0.5 .* (acceleration .+ new_accel) .* dtime
end

function get_axlen(all_bodies::Vector{Body})
    central_pos = all_bodies[1].position
    max_dist = 0.0
    for body in all_bodies[2:end]
        dist = norm(body.position .- central_pos)
        max_dist = max(max_dist, dist)
    end
    return max_dist*1.1
end


function simulate(all_bodies::Vector{Body}, time::Float64, steps::Int64, frames::Int64, name::String)
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

function simulate_without_save(all_bodies::Vector{Body}, time::Float64, steps::Int64)
    for _ in 1:steps
        for body in all_bodies
            update!(body, all_bodies, time / steps)
        end
    end
end

end # module BasicNBodySim
