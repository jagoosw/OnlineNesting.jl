# An M2 tide on the shelf, approaches and bay grids, nested one way: each nest takes its
# boundary state from its parent and nothing is fed back. The tide is a uniform sea-level
# oscillation on the four open sides of the shelf, not a tidal atlas. The grids are twice as
# coarse as the target case (5 km, 1 km and 200 m) with four layers, to run in minutes.
#
#   julia --project -t auto validation/pembrokeshire_tide.jl [M2 periods, default 2]
#
# For a GPU, load CUDA before this script and set `arch = GPU()`.

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: TracerReservoir
using Oceananigans.Grids: Open, static_column_depthᶜᶠᵃ, static_column_depthᶠᶜᵃ, λnodes, φnodes
using Oceananigans.Operators: ℑxyᶜᶠᵃ, ℑxyᶠᶜᵃ
using NumericalEarth
using NumericalEarth.DataWrangling: BoundingBox, Metadatum
using OnlineNesting
using CairoMakie
using Printf

arch = CPU()
periods = isempty(ARGS) ? 2 : parse(Float64, ARGS[1])
tidal_period = 12.4206hours
Δt = 200 # s, on the shelf
refinement = 5
width = 3 # parent cells
minimum_depth = 10
z = [-250, -150, -95, -40, 0]

#####
##### Grids and bathymetry
#####

shelf_grid = LatitudeLongitudeGrid(arch; size=(110, 122, 4), longitude=(-10.5, -2.5), latitude=(49, 54.5), halo=(7, 7, 4),
                                   topology=(Open, Open, Bounded), z)
approaches_grid = nested_grid(shelf_grid; longitude=(-6.4, -4.6), latitude=(51.3, 52.4), depth=150, refinement)
bay_grid = nested_grid(approaches_grid; longitude=(-5.50, -5.06), latitude=(51.68, 51.94), depth=95, refinement)

function emodnet_bottom_height(grid; interpolation_passes, major_basins=Inf)
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(grid))
    bottom_height = regrid_bathymetry(grid, metadatum; interpolation_passes, major_basins)
    return refer_to_mean_sea_level!(bottom_height, 0; minimum_depth)
end

shelf_bottom = emodnet_bottom_height(shelf_grid; interpolation_passes=4, major_basins=1)
approaches_bottom = emodnet_bottom_height(approaches_grid; interpolation_passes=2)
bay_bottom = emodnet_bottom_height(bay_grid; interpolation_passes=1)

blend_parent_bottom!(approaches_bottom, shelf_bottom; width, minimum_depth)
blend_parent_bottom!(bay_bottom, approaches_bottom; width, minimum_depth)
restrict_bottom!(approaches_bottom, bay_bottom; width, minimum_depth)
restrict_bottom!(shelf_bottom, approaches_bottom; width, minimum_depth)

sea = map(bottom -> Array(interior(bottom, :, :, 1)) .< 0, (shelf = shelf_bottom, approaches = approaches_bottom, bay = bay_bottom))

immersed(grid, bottom) = ImmersedBoundaryGrid(grid, PartialCellBottom(Array(interior(bottom, :, :, 1))))
shelf_grid, approaches_grid, bay_grid = immersed.((shelf_grid, approaches_grid, bay_grid), (shelf_bottom, approaches_bottom, bay_bottom))

#####
##### Models
#####

initial_temperature(z) = 10 + 8 * exp(z / 50)

# Quadratic drag on the depth-mean flow: -Cᴰ |u| u / H
@inline function u_drag(i, j, k, grid, clock, fields, Cᴰ)
    u = @inbounds fields.u[i, j, k]
    v = ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    H = static_column_depthᶠᶜᵃ(i, j, grid)
    return ifelse(H > 0, -Cᴰ * sqrt(u^2 + v^2) * u / H, zero(grid))
end

@inline function v_drag(i, j, k, grid, clock, fields, Cᴰ)
    u = ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)
    v = @inbounds fields.v[i, j, k]
    H = static_column_depthᶜᶠᵃ(i, j, grid)
    return ifelse(H > 0, -Cᴰ * sqrt(u^2 + v^2) * v / H, zero(grid))
end

drag = (u = Forcing(u_drag; discrete_form=true, parameters=2.5e-3), v = Forcing(v_drag; discrete_form=true, parameters=2.5e-3))

function ocean(grid; boundary_conditions, forcing = drag)
    spacing = first(diff(Array(φnodes(grid, Face())))) * π * grid.radius / 180
    closure = (HorizontalScalarDiffusivity(ν = 0.01 * spacing),
               VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=1e-3, κ=1e-5))

    model = HydrostaticFreeSurfaceModel(grid; boundary_conditions, forcing, closure,
                                        tracers = :T,
                                        buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(), constant_salinity=35),
                                        coriolis = HydrostaticSphericalCoriolis(),
                                        momentum_advection = WENOVectorInvariant(),
                                        tracer_advection = WENO(),
                                        free_surface = SplitExplicitFreeSurface(grid; substeps=20),
                                        timestepper = :SplitRungeKutta3)

    set!(model, T = (λ, φ, z) -> initial_temperature(z))
    return model
end

function nested_ocean(parent, grid)
    exchanger = StateExchanger(parent)
    sponge = nested_forcings(exchanger, grid; rate=1/1hour, width=refinement*width)
    forcing = (u = (sponge.u, drag.u), v = (sponge.v, drag.v), T = sponge.T)
    return ocean(grid; boundary_conditions=nested_boundary_conditions(exchanger, grid), forcing), exchanger
end

struct Tide{FT}
    amplitude :: FT
    frequency :: FT
    ramp :: FT
end

@inline (tide::Tide)(i, j, grid, clock, fields) = tide.amplitude * sin(tide.frequency * clock.time) * tanh(clock.time / tide.ramp)

function nested_levels(; amplitude)
    radiation = GravityWaveRadiationBoundaryCondition((0, Tide(amplitude, 2π / tidal_period, 3hours)))
    wall = NormalFlowBoundaryCondition(0)
    reservoir = ValueBoundaryCondition((ξ, z, t) -> initial_temperature(z); scheme=TracerReservoir())

    boundary_conditions = (U = FieldBoundaryConditions(west=radiation, east=radiation),
                           V = FieldBoundaryConditions(south=radiation, north=radiation),
                           u = FieldBoundaryConditions(west=wall, east=wall),
                           v = FieldBoundaryConditions(south=wall, north=wall),
                           T = FieldBoundaryConditions(west=reservoir, east=reservoir, south=reservoir, north=reservoir))

    shelf = ocean(shelf_grid; boundary_conditions)
    approaches, approaches_exchanger = nested_ocean(shelf, approaches_grid)
    bay, bay_exchanger = nested_ocean(approaches, bay_grid)

    inner_nest = TwoWayNestedModel(approaches, bay, bay_exchanger; refinement)
    return TwoWayNestedModel(shelf, inner_nest, approaches_exchanger; refinement), (; shelf, approaches, bay)
end

#####
##### Diagnostics
#####

sea_level(model) = Array(interior(model.free_surface.displacement, :, :, 1))

function surface_speed(model)
    Nz = size(model.grid, 3)
    u = Array(interior(model.velocities.u, :, :, Nz))
    v = Array(interior(model.velocities.v, :, :, Nz))
    return sqrt.(((u[1:end-1, :] .+ u[2:end, :]) ./ 2) .^ 2 .+ ((v[:, 1:end-1] .+ v[:, 2:end]) ./ 2) .^ 2)
end

# The column of `grid` containing (λ, φ), or `nothing` if it is outside the grid or on land
function column(grid, sea, (λ, φ))
    λᶠ, φᶠ = Array(λnodes(grid, Face())), Array(φnodes(grid, Face()))
    (λᶠ[1] ≤ λ < λᶠ[end] && φᶠ[1] ≤ φ < φᶠ[end]) || return nothing
    i, j = searchsortedlast(λᶠ, λ), searchsortedlast(φᶠ, φ)
    return sea[i, j] ? (i, j) : nothing
end

probe(values, ::Nothing) = NaN
probe(values, column) = values[column...]

#####
##### Control: nothing imposed
#####

still_nest, still = nested_levels(amplitude = 0.0)
for _ in 1:3
    time_step!(still_nest, Δt)
end

println("\nControl, no tide imposed, 3 shelf steps:")
for (name, model) in pairs(still)
    @printf("  %-12s max|u| = %.2e, max|v| = %.2e m s⁻¹, max|η| = %.2e m\n", name,
            maximum(abs, model.velocities.u), maximum(abs, model.velocities.v), maximum(abs, model.free_surface.displacement))
end

#####
##### The tide
#####

nest, levels = nested_levels(amplitude = 1.0)
places = (st_brides_bay = (-5.25, 51.80), ramsey_sound = (-5.3225, 51.868))
columns = map((model, sea) -> map(place -> column(model.grid, sea, place), places), levels, sea)
series = (time = Float64[], η = map(_ -> map(_ -> Float64[], places), levels), speed = map(_ -> map(_ -> Float64[], places), levels))

println("\nLevels:")
for (n, (name, model)) in enumerate(pairs(levels))
    @printf("  %-12s %3d × %3d × %d, %2.0f%% sea, Δt = %g s\n", name, size(model.grid)..., 100 * count(sea[name]) / length(sea[name]),
            Δt / refinement^(n - 1))
end

# Surface speed and sea level on every level, with land blanked
function surface_fields(name)
    speed, η = surface_speed(levels[name]), sea_level(levels[name])
    speed[.!sea[name]] .= NaN
    η[.!sea[name]] .= NaN
    return speed, η
end

figure = Figure(size=(1500, 1200))
title = Observable("")
Label(figure[0, 1:3], title; halign=:left, fontsize=17)

panels = map(enumerate(keys(levels))) do (n, name)
    grid = levels[name].grid
    λ, φ = Array(λnodes(grid, Center())), Array(φnodes(grid, Center()))
    aspect = AxisAspect((λ[end] - λ[1]) * cosd((φ[1] + φ[end]) / 2) / (φ[end] - φ[1]))
    speed, η = map(Observable, surface_fields(name))
    speed_limit = (1.5, 1.5, 2.5)[n]

    axis = Axis(figure[1, n]; title = @sprintf("%s, %d × %d × %d", name, size(grid)...), aspect, ylabel = n == 1 ? "latitude (°N)" : "")
    heatmap!(axis, λ, φ, speed; colormap=:Blues, colorrange=(0, speed_limit), nan_color=:gray80)
    Colorbar(figure[2, n]; colormap=:Blues, limits=(0, speed_limit), vertical=false, flipaxis=false, label="surface speed (m s⁻¹)")

    axis = Axis(figure[3, n]; aspect, xlabel="longitude (°E)", ylabel = n == 1 ? "latitude (°N)" : "")
    heatmap!(axis, λ, φ, η; colormap=:balance, colorrange=(-3, 3), nan_color=:gray80)
    Colorbar(figure[4, n]; colormap=:balance, limits=(-3, 3), vertical=false, flipaxis=false, label="sea level (m)")

    (; speed, η)
end

function update_panels!()
    for (panel, name) in zip(panels, keys(levels))
        panel.speed[], panel.η[] = surface_fields(name)
    end
    title[] = @sprintf("One-way nested M2 tide after %.1f h", nest.clock.time / 3600)
end

movie = VideoStream(figure; framerate=15)
steps = round(Int, periods * tidal_period / Δt)
start = time()
println()

for step in 1:steps
    time_step!(nest, Δt)
    push!(series.time, nest.clock.time)

    for name in keys(levels), place in keys(places)
        push!(series.η[name][place], probe(sea_level(levels[name]), columns[name][place]))
        push!(series.speed[name][place], probe(surface_speed(levels[name]), columns[name][place]))
    end

    if step % 3 == 0
        update_panels!()
        recordframe!(movie)
    end

    if step % max(1, steps ÷ 12) == 0 || step == steps
        @printf("  t = %5.1f h  |  %s  |  wall time %.0f s\n", nest.clock.time / 3600,
                join([@sprintf("%s: max|η| %.2f m, max speed %.2f m/s", name, maximum(abs, sea_level(model)), maximum(surface_speed(model)))
                      for (name, model) in pairs(levels)], "  |  "), time() - start)
        any(model -> any(isnan, sea_level(model)), levels) && error("NaN in the sea level")
    end
end

#####
##### Agreement of each nest with its parent
#####

late = series.time .> tidal_period * min(1, periods / 2)
rms(a) = sqrt(sum(abs2, a) / length(a))
correlation(a, b) = sum((a .- sum(a) / length(a)) .* (b .- sum(b) / length(b))) / (rms(a .- sum(a) / length(a)) * rms(b .- sum(b) / length(b)) * length(a))

println("\nEach nest against its parent at the same place, after the first ", round(tidal_period * min(1, periods / 2) / 3600, digits=1), " h:")
for (parent, child) in ((:shelf, :approaches), (:approaches, :bay)), place in keys(places)
    (isnothing(columns[parent][place]) || isnothing(columns[child][place])) && continue
    η, parent_η = series.η[child][place][late], series.η[parent][place][late]
    speed, parent_speed = series.speed[child][place][late], series.speed[parent][place][late]
    @printf("  %-10s vs %-10s at %-13s sea level amplitude %.2f vs %.2f m, rms difference %.3f m, correlation %.4f\n",
            child, parent, place, maximum(abs, η), maximum(abs, parent_η), rms(η .- parent_η), correlation(η, parent_η))
    @printf("  %-10s vs %-10s at %-13s surface speed peak %.2f vs %.2f m/s, rms difference %.3f m/s\n",
            child, parent, place, maximum(speed), maximum(parent_speed), rms(speed .- parent_speed))
end

#####
##### Movie and figure
#####

save(joinpath(@__DIR__, "pembrokeshire_tide.mp4"), movie)
update_panels!()

hours = series.time ./ 3600

sea_level_axis = Axis(figure[5, 1:3]; ylabel="sea level, St Brides Bay (m)")
speed_axis = Axis(figure[6, 1:3]; xlabel="time (h)", ylabel="surface speed, Ramsey Sound (m s⁻¹)")

for name in keys(levels)
    lines!(sea_level_axis, hours, series.η[name].st_brides_bay; label=string(name))
    lines!(speed_axis, hours, series.speed[name].ramsey_sound; label=string(name))
end

axislegend(sea_level_axis; position=:lt)
save(joinpath(@__DIR__, "pembrokeshire_tide.png"), figure)
