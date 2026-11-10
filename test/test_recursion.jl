using Test
using Oceananigans
using Oceananigans.Units
using Oceananigans.TimeSteppers: update_state!
using NumericalEarth
using OnlineNesting

function ocean(grid; kwargs...)
    return HydrostaticFreeSurfaceModel(grid; free_surface = SplitExplicitFreeSurface(grid; substeps=20),
                                       timestepper = :SplitRungeKutta3,
                                       tracers = (:T, :S, :c),
                                       buoyancy = SeawaterBuoyancy(),
                                       kwargs...)
end

function nested_ocean(parent, grid)
    exchanger = StateExchanger(parent)
    child = ocean(grid; boundary_conditions = nested_boundary_conditions(exchanger, grid),
                        forcing = nested_forcings(exchanger, grid; rate=1/1hour, width=4))
    return child, exchanger
end

function three_levels()
    shelf_grid = LatitudeLongitudeGrid(size=(40, 40, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0))
    approaches_grid = nested_grid(shelf_grid; longitude=(3, 7), latitude=(43, 47), refinement=3)
    bay_grid = nested_grid(approaches_grid; longitude=(4, 6), latitude=(44, 46), refinement=3)

    shelf = ocean(shelf_grid)
    approaches, approaches_exchanger = nested_ocean(shelf, approaches_grid)
    bay, bay_exchanger = nested_ocean(approaches, bay_grid)

    inner_nest = TwoWayNestedModel(approaches, bay, bay_exchanger; width=2)
    return TwoWayNestedModel(shelf, inner_nest, approaches_exchanger; width=2), (shelf, approaches, bay)
end

stratification(λ, φ, z) = 10 + z / 20
bump(λ, φ, z) = 0.1 * exp(-((λ - 2)^2 + (φ - 45)^2) / 2)

@testset "A gravity wave crosses three levels" begin
    nest, levels = three_levels()
    shelf, approaches, bay = levels
    foreach(level -> set!(level, T=stratification, S=35, c=1, η=bump), levels)

    i, j, _ = parent_indices(approaches.grid, bay.grid)
    block = (3, length(i), 3, length(j))
    amplitude, difference = 0.0, 0.0
    synchronized = true

    for _ in 1:80
        time_step!(nest, 5minutes)
        synchronized &= shelf.clock.time == approaches.clock.time == bay.clock.time

        η = Array(interior(approaches.free_surface.displacement))[i, j, 1]
        blocks = reshape(Array(interior(bay.free_surface.displacement)), block)
        restricted_η = sum(blocks, dims=(1, 3))[1, :, 1, :] ./ 9
        amplitude = max(amplitude, maximum(abs, η))
        difference = max(difference, maximum(abs, restricted_η .- η))
    end

    @test synchronized
    @test bay.clock.iteration == 3 * approaches.clock.iteration == 9 * shelf.clock.iteration
    @test difference < 0.05 * amplitude

    for level in levels
        @test maximum(abs, Array(interior(level.tracers.c)) .- 1) < 1e-12
    end
end

@testset "The clocks of three levels agree with a step that is not representable" begin
    nest, levels = three_levels()
    foreach(level -> set!(level, T=stratification, S=35), levels)
    synchronized = true

    for _ in 1:20
        time_step!(nest, 0.1)
        synchronized &= allequal(level.clock.time for level in levels)
    end

    @test synchronized
end

@testset "The innermost level reaches the outermost at initialization" begin
    nest, levels = three_levels()
    shelf, approaches, bay = levels
    foreach(level -> set!(level, T=stratification, S=35), levels)
    set!(bay, c=1)
    update_state!(nest)

    @test maximum(interior(approaches.tracers.c)) ≈ 1
    @test maximum(interior(shelf.tracers.c)) ≈ 1
end

function earth_system_ocean(grid; kwargs...)
    ocean = ocean_simulation(grid; free_surface = SplitExplicitFreeSurface(grid; substeps=20),
                             closure = nothing,
                             coriolis = nothing,
                             momentum_advection = VectorInvariant(),
                             tracer_advection = Centered(),
                             kwargs...)
    return OceanOnlyModel(ocean)
end

@testset "Earth system models as levels stay at rest" begin
    parent_grid = LatitudeLongitudeGrid(size=(40, 40, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0))
    grid = nested_grid(parent_grid; longitude=(3, 7), latitude=(43, 47), refinement=3)

    parent = earth_system_ocean(parent_grid)
    exchanger = StateExchanger(parent)
    child = earth_system_ocean(grid; boundary_conditions = nested_boundary_conditions(exchanger, grid),
                                     forcing = nested_forcings(exchanger, grid; rate=1/1hour, width=4))

    oceans = (parent.ocean.model, child.ocean.model)
    foreach(ocean -> set!(ocean, T=stratification, S=35), oceans)

    simulation = Simulation(TwoWayNestedModel(parent, child, exchanger; width=2); Δt=5minutes, stop_iteration=6, verbose=false)
    run!(simulation)

    @test parent.clock.time == child.clock.time == 30minutes
    @test all(ocean.clock.time == 30minutes for ocean in oceans)
    @test child.ocean.model.clock.iteration == 3 * parent.ocean.model.clock.iteration

    for ocean in oceans
        @test maximum(abs, ocean.velocities.u) < 1e-12
        @test maximum(abs, ocean.free_surface.displacement) < 1e-12
    end
end

@testset "Feedback refreshes the ocean of an Earth system model" begin
    parent_grid = LatitudeLongitudeGrid(size=(20, 20, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0))
    grid = nested_grid(parent_grid; longitude=(3, 7), latitude=(43, 47), refinement=3)

    parent = earth_system_ocean(parent_grid)
    exchanger = StateExchanger(parent)
    child = earth_system_ocean(grid; boundary_conditions = nested_boundary_conditions(exchanger, grid))

    set!(parent.ocean.model, T=stratification, S=35)
    set!(child.ocean.model, T=(λ, φ, z) -> stratification(λ, φ, z) + 1, S=35)
    update_state!(TwoWayNestedModel(parent, child, exchanger; width=2))

    pressure = Array(interior(parent.ocean.model.pressure.pHY′))
    update_state!(parent.ocean.model)
    @test pressure == Array(interior(parent.ocean.model.pressure.pHY′))
    @test maximum(abs, pressure) > 0
end
