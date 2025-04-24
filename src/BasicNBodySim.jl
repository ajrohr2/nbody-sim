module BasicNBodySim

using LinearAlgebra

export simulate, update!, get_axlen, Body, Planet, Moon, Star

abstract type Body end

"""
    Planet(mass::Float64, velocity::Vector{Float64}, position::Vector{Float64}, color) <: Body

Create a Planet object with `mass`, `velocity`, `position`, and `color` in S.I. units.

# Examples
```jldoctest; output = false
julia> Planet(5.9722e24, [0.0, 0.0], [0.0, 0.0], :blue)
Planet(5.9722e24, [0.0, 0.0], [0.0, 0.0], :blue)
```
"""
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


"""
    Star(mass::Float64, velocity::Vector{Float64}, position::Vector{Float64}, color) <: Body

Create a Star object with `mass`, `velocity`, `position`, and `color` in S.I. units.

# Examples
```jldoctest; output = false
julia> Star(1.989e30, [0.0, 0.0], [0.0, 0.0], :blue)
Star(1.989e30, [0.0, 0.0], [0.0, 0.0], :blue)
```
"""
mutable struct Star <: Body
    mass:: Float64
    velocity::Vector{Float64}
    position::Vector{Float64}
    color::Symbol
    function Star(mass::Float64, velocity::Vector{Float64}, position::Vector{Float64}, color)
        new(min(mass, 3e32), velocity, position, color)
    end
end


"""
    Moon(mass::Float64, velocity::Vector{Float64}, position::Vector{Float64}, color) <: Body

Create a Moon object with `mass`, `velocity`, `position`, and `color` in S.I. units.

# Examples
```jldoctest; output = false
julia> Moon(7.348e22, [3.844e8, 0.0], [0.0, 1022.0], :gray)
Moon(7.348e22, [3.844e8, 0.0], [0.0, 1022.0], :gray)
```
"""
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

"""
    simulate(all_bodies::Vector{Body}, time::Float64, steps::Int64)

Simulate the gravitational force between bodies in `all_bodies` for `time` seconds using Verlet Velocity symplectic integration.

The delta time is determined by `time / steps`.
"""
function simulate(all_bodies::Vector{Body}, time::Float64, steps::Int64)
    for _ in 1:steps
        for body in all_bodies
            update!(body, all_bodies, time / steps)
        end
    end
end

end # module BasicNBodySim
