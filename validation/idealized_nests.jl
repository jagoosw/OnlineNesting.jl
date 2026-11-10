# Videos of two idealized three-level nests, shelf → approaches → bay, each with flat 100 m
# bottoms and refinement 3:
# - a tracer blob carried around a periodic channel, through the nests and back into them;
# - a gravity wave from a sea-level bump in a closed box, crossing the nests;
# - a geostrophic eddy set in the bay only, one shelf cell in radius.
# The left panel is the shelf alone, with the state fed back from the nests; the right panel
# draws each nest over its parent.
#
#   julia --project -t auto validation/idealized_nests.jl

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: x_domain, y_domain, λnodes, φnodes
using OnlineNesting
using CairoMakie

function ocean(grid; kwargs...)
    return HydrostaticFreeSurfaceModel(grid; free_surface = SplitExplicitFreeSurface(grid; substeps=20),
                                       timestepper = :SplitRungeKutta3,
                                       tracers = (:T, :S, :c),
                                       buoyancy = SeawaterBuoyancy(),
                                       kwargs...)
end

function nested_ocean(parent, grid; kwargs...)
    exchanger = StateExchanger(parent)
    child = ocean(grid; boundary_conditions = nested_boundary_conditions(exchanger, grid),
                        forcing = nested_forcings(exchanger, grid; rate=1/1hour, width=4), kwargs...)
    return child, exchanger
end

function three_levels(shelf_grid; longitude, latitude, kwargs...)
    approaches_grid = nested_grid(shelf_grid; longitude=longitude[1], latitude=latitude[1], refinement=3)
    bay_grid = nested_grid(approaches_grid; longitude=longitude[2], latitude=latitude[2], refinement=3)

    shelf = ocean(shelf_grid; kwargs...)
    approaches, approaches_exchanger = nested_ocean(shelf, approaches_grid; kwargs...)
    bay, bay_exchanger = nested_ocean(approaches, bay_grid; kwargs...)

    inner_nest = TwoWayNestedModel(approaches, bay, bay_exchanger; width=2)
    return TwoWayNestedModel(shelf, inner_nest, approaches_exchanger; width=2), (shelf, approaches, bay)
end

function outline!(axis, grid)
    (west, east), (south, north) = x_domain(grid), y_domain(grid)
    return lines!(axis, [west, east, east, west, west], [south, south, north, north, south]; color=:black, linewidth=1)
end

# `surface_field(level)` is the horizontal array to draw
function record_levels(surface_field, nest, levels, filename; Δt, steps, title, colormap, colorrange)
    shelf_grid = first(levels).grid
    aspect = DataAspect()
    figure = Figure(size=(1100, 520))
    time_label = Observable("")
    Label(figure[0, 1:2], time_label)
    shelf_axis = Axis(figure[1, 1]; title="shelf", xlabel="longitude (°E)", ylabel="latitude (°N)", aspect)
    composite_axis = Axis(figure[1, 2]; title="shelf, approaches and bay", xlabel="longitude (°E)", aspect)

    fields = [Observable(Array(surface_field(level))) for level in levels]
    nodes(level) = (λnodes(level.grid, Center()), φnodes(level.grid, Center()))

    heatmap!(shelf_axis, nodes(levels[1])..., fields[1]; colormap, colorrange)
    for (level, field) in zip(levels, fields)
        heatmap!(composite_axis, nodes(level)..., field; colormap, colorrange)
    end
    for level in levels[2:end], axis in (shelf_axis, composite_axis)
        outline!(axis, level.grid)
    end
    Colorbar(figure[1, 3]; colormap, colorrange, label=title)

    record(figure, filename, 0:steps; framerate=24) do n
        n > 0 && time_step!(nest, Δt)
        for (level, field) in zip(levels, fields)
            field[] = Array(surface_field(level))
        end
        time_label[] = "$title after $(prettytime(first(levels).clock.time))"
    end

    return filename
end

stratification(λ, φ, z) = 10 + z / 20

#####
##### Blob
#####

channel_grid = LatitudeLongitudeGrid(size=(48, 16, 4), longitude=(0, 12), latitude=(43, 47), z=(-100, 0),
                                     topology=(Periodic, Bounded, Bounded))
nest, levels = three_levels(channel_grid; longitude=((4, 8), (5, 7)), latitude=((44, 46), (44.5, 45.5)))

blob(λ, φ, z) = exp(-((λ - 6)^2 + (φ - 45)^2) / (2 * 0.3^2))
flow(λ, φ, z) = 4 * cosd(φ)
foreach(level -> set!(level, T=stratification, S=35, c=blob, u=flow), levels)

Δt = 5minutes
circuit = round(Int, deg2rad(12) * channel_grid.radius / 4 / Δt)
record_levels(level -> interior(level.tracers.c, :, :, size(level.grid, 3)), nest, levels,
              joinpath(@__DIR__, "idealized_blob.mp4");
              Δt = 2Δt, steps = circuit ÷ 2, title = "surface tracer", colormap = :dense, colorrange = (0, 1))

#####
##### Wave
#####

box_grid = LatitudeLongitudeGrid(size=(40, 40, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0))
nest, levels = three_levels(box_grid; longitude=((3, 7), (4, 6)), latitude=((43, 47), (44, 46)))

bump(λ, φ, z) = 0.1 * exp(-((λ - 2)^2 + (φ - 45)^2) / 2)
foreach(level -> set!(level, T=stratification, S=35, η=bump), levels)

record_levels(level -> interior(level.free_surface.displacement, :, :, 1), nest, levels,
              joinpath(@__DIR__, "idealized_wave.mp4");
              Δt, steps = 160, title = "sea level (m)", colormap = :balance, colorrange = (-0.03, 0.03))

#####
##### Eddy
#####

nest, levels = three_levels(box_grid; longitude=((3, 7), (4, 6)), latitude=((43, 47), (44, 46)),
                            coriolis = HydrostaticSphericalCoriolis())
shelf, approaches, bay = levels

radius = 0.25
g = shelf.free_surface.gravitational_acceleration
fR(φ) = 2 * shelf.coriolis.rotation_rate * sind(φ) * deg2rad(box_grid.radius)
η(λ, φ, z=0) = 0.02 * exp(-(((λ - 5) * cosd(45))^2 + (φ - 45)^2) / 2radius^2)
u(λ, φ, z) =  g / fR(φ) * η(λ, φ) * (φ - 45) / radius^2
v(λ, φ, z) = -g / fR(φ) * η(λ, φ) * (λ - 5) * cosd(45)^2 / (cosd(φ) * radius^2)

foreach(level -> set!(level, T=stratification, S=35), levels)
set!(bay; η, u, v)

record_levels(level -> interior(level.free_surface.displacement, :, :, 1), nest, levels,
              joinpath(@__DIR__, "idealized_eddy.mp4");
              Δt, steps = 144, title = "sea level (m)", colormap = :balance, colorrange = (-0.02, 0.02))
