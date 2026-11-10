using Test
using Oceananigans
using Oceananigans.Units
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Operators: Azᶜᶜᶜ
using Oceananigans.TimeSteppers: update_state!
using NumericalEarth
using OnlineNesting

function ocean(grid; kwargs...)
    return HydrostaticFreeSurfaceModel(grid; free_surface = SplitExplicitFreeSurface(grid; substeps=20),
                                       timestepper = :SplitRungeKutta3,
                                       tracers = (:T, :S),
                                       buoyancy = SeawaterBuoyancy(),
                                       kwargs...)
end

diffusive_ocean(grid; kwargs...) = ocean(grid; tracer_advection = WENO(order=5), closure = HorizontalScalarDiffusivity(κ=500), kwargs...)

function earth_system_ocean(grid; kwargs...)
    ocean = ocean_simulation(grid; free_surface = SplitExplicitFreeSurface(grid; substeps=20),
                             closure = nothing,
                             coriolis = nothing,
                             momentum_advection = VectorInvariant(),
                             tracer_advection = Centered(),
                             kwargs...)
    return OceanOnlyModel(ocean)
end

refinement, width = 3, 2

# A channel with a nest over each of `extents`, every nest inside the one before, on the grids `immerse`
# returns for the list of level grids: the nested model and its levels
function nested_channel(level, z, extents; immerse = identity)
    grids = Any[LatitudeLongitudeGrid(size=(24, 16, 4), longitude=(0, 6), latitude=(43, 47), halo=(5, 5, 5),
                                      topology=(Periodic, Bounded, Bounded); z)]

    for extent in extents
        push!(grids, nested_grid(grids[end]; extent..., refinement))
    end

    grids = immerse(grids)
    levels = Any[level(grids[1])]
    exchangers = []

    for grid in grids[2:end]
        exchanger = StateExchanger(levels[end])
        push!(exchangers, exchanger)
        push!(levels, level(grid; boundary_conditions = nested_boundary_conditions(exchanger, grid),
                                  forcing = nested_forcings(exchanger, grid; rate=1/1hour, width=4)))
    end

    nest = foldr(zip(levels[1:end-1], exchangers); init=levels[end]) do (parent, exchanger), child
        TwoWayNestedModel(parent, child, exchanger; width)
    end

    return nest, levels
end

cell_measure(measure, grid) = Array(interior(compute!(Field(KernelFunctionOperation{Center, Center, Center}(measure, grid)))))
heat(model) = cell_measure(OnlineNesting.wet_volume, model.grid) .* Array(interior(model.tracers.T))
volume(model) = cell_measure(Azᶜᶜᶜ, model.grid)[:, :, 1] .* Array(interior(model.free_surface.displacement))[:, :, 1]

# The sum of `measure` over the cells of each of `models`, finest last, that no nest feeds back to
function composite(measure, models)
    i, j = Colon(), Colon()
    total = 0.0

    for (parent, child) in zip(models[1:end-1], models[2:end])
        values = measure(parent)
        nest_i, nest_j, _ = parent_indices(parent.grid, child.grid)
        total += sum(values[i, j, :]) - sum(values[nest_i[1+width:end-width], nest_j[1+width:end-width], :])

        Nx, Ny, _ = size(child.grid)
        rim = refinement * width
        i, j = 1+rim:Nx-rim, 1+rim:Ny-rim
    end

    return total + sum(measure(models[end])[i, j, :])
end

# a warm patch and a sea level bump in a sheared flow around the channel
temperature(λ, φ, z) = 10 + z / 20 + exp(-((λ - 3)^2 + (φ - 45)^2) / (2 * 0.3^2)) * (1 - z / 100)
flow(λ, φ, z) = 4 * cosd(φ) * (1 + z / 100)
bump(λ, φ, z) = 0.1 * exp(-((λ - 2)^2 + (φ - 45)^2) / 2)

salinity_departure(model) = maximum(abs, Array(interior(model.tracers.S))[cell_measure(OnlineNesting.wet_volume, model.grid) .> 0] .- 35)

# The largest changes of the composite heat and volume, relative to the heat and to the volume
# of the bump, and the largest departure of the salinity from uniform
function conservation_errors(nest, models; steps, spin_up = Returns(nothing))
    foreach(model -> set!(model, T=temperature, S=35, u=flow, η=bump), models)
    spin_up()
    update_state!(nest)

    initial_heat = composite(heat, models)
    initial_volume = composite(volume, models)
    bump_volume = sum(abs, volume(models[1]))
    heat_error, volume_error, salinity_error = 0.0, 0.0, 0.0

    for _ in 1:steps
        time_step!(nest, 5minutes)
        heat_error = max(heat_error, abs(composite(heat, models) / initial_heat - 1))
        volume_error = max(volume_error, abs(composite(volume, models) - initial_volume) / bump_volume)
        salinity_error = max(salinity_error, maximum(salinity_departure, models))
    end

    return heat_error, volume_error, salinity_error
end

z_star() = MutableVerticalDiscretization((-100, 0))
approaches = (longitude=(2, 4), latitude=(44, 46))
bay = (longitude=(2.5, 3.5), latitude=(44.5, 45.5))

# a shelf over the western rim of the approaches and the parent cells beside its feedback region
shelf(λ, φ) = λ < 2.5 ? -50 : -100
parent_shelf(grids) = [ImmersedBoundaryGrid(grids[1], GridFittedBottom(shelf)), grids[2]]
child_shelf(grids) = [grids[1], ImmersedBoundaryGrid(grids[2], GridFittedBottom(shelf))]

static_cases = (("a flat bottom", identity), ("a parent shelf beside a deeper child", parent_shelf))

@testset "Two levels on a static grid with $name conserve volume and keep a uniform tracer" for (name, immerse) in static_cases
    nest, models = nested_channel(ocean, (-100, 0), (approaches,); immerse)
    _, volume_error, salinity_error = conservation_errors(nest, models; steps=40)

    @test volume_error < 1e-12
    @test salinity_error < 1e-11
end

# a sloping sea floor, and on the nest bumps the parent does not resolve
sloping_bottom(λ, φ) = -30 - 10λ
bumpy_bottom(λ, φ) = sloping_bottom(λ, φ) + 8 * sinpi(4λ) * sinpi(4φ)

function reconciled_bottoms((parent_grid, grid))
    parent_bottom_height = Field{Center, Center, Nothing}(parent_grid)
    bottom_height = Field{Center, Center, Nothing}(grid)
    set!(parent_bottom_height, sloping_bottom)
    set!(bottom_height, bumpy_bottom)

    blend_parent_bottom!(bottom_height, parent_bottom_height; width, minimum_depth=10)
    restrict_bottom!(parent_bottom_height, bottom_height; width, minimum_depth=10)

    return [ImmersedBoundaryGrid(parent_grid, PartialCellBottom(parent_bottom_height)),
            ImmersedBoundaryGrid(grid, PartialCellBottom(bottom_height))]
end

z_star_cases = (("two levels with upwind advection and diffusion", diffusive_ocean, identity, (approaches,), identity, 40),
                ("three levels", ocean, identity, (approaches, bay), identity, 20),
                ("Earth system models", earth_system_ocean, level -> level.ocean.model, (approaches,), identity, 10),
                ("two levels with a child shelf beside a deeper parent", ocean, identity, (approaches,), child_shelf, 20),
                ("two levels with a parent shelf beside a deeper child", ocean, identity, (approaches,), parent_shelf, 20),
                ("two levels with reconciled bottoms", ocean, identity, (approaches,), reconciled_bottoms, 20))

@testset "On z-star grids, $name conserve volume and tracer" for (name, level, ocean_model, extents, immerse, steps) in z_star_cases
    nest, levels = nested_channel(level, z_star(), extents; immerse)
    heat_error, volume_error, salinity_error = conservation_errors(nest, ocean_model.(levels); steps)

    @test heat_error < 1e-12
    @test volume_error < 1e-12
    @test salinity_error < 1e-11
end

@testset "A parent stepped alone leaves the next steps of the nest conservative" begin
    nest, levels = nested_channel(earth_system_ocean, z_star(), (approaches,))
    models = [level.ocean.model for level in levels]
    heat_error, volume_error, _ = conservation_errors(nest, models; steps=5, spin_up = () -> time_step!(levels[1], 5minutes))

    @test heat_error < 1e-12
    @test volume_error < 1e-12
end
