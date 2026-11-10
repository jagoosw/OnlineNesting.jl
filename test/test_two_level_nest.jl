using Test
using Oceananigans
using Oceananigans.Units
using OnlineNesting

function ocean(grid; kwargs...)
    return HydrostaticFreeSurfaceModel(grid; free_surface = SplitExplicitFreeSurface(grid; substeps=20),
                                       timestepper = :SplitRungeKutta3,
                                       tracers = (:T, :S),
                                       buoyancy = SeawaterBuoyancy(),
                                       kwargs...)
end

function nested_levels(; refinement, parent_bottom = -100)
    underlying_grid = LatitudeLongitudeGrid(size=(40, 40, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0))
    parent_grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(parent_bottom))
    grid = nested_grid(parent_grid; longitude=(3, 7), latitude=(43, 47), refinement)

    parent = ocean(parent_grid)
    exchanger = StateExchanger(parent)
    child = ocean(grid; boundary_conditions = nested_boundary_conditions(exchanger, grid),
                        forcing = nested_forcings(exchanger, grid; rate=1/1hour, width=4))

    return parent, child, TwoWayNestedModel(parent, child, exchanger; width=2)
end

stratification(λ, φ, z) = 10 + z / 20

# land along the west side of the nest and a shelf 50 m deep along its south side
coast(λ, φ) = ifelse(λ < 3, 0, ifelse(φ < 43.5, -50, -100))

@testset "A stratified ocean at rest stays at rest, parent bottom $parent_bottom" for parent_bottom in (-100, coast)
    parent, child, nest = nested_levels(; refinement=3, parent_bottom)
    set!(parent, T=stratification, S=35)
    set!(child, T=stratification, S=35)

    simulation = Simulation(nest; Δt=5minutes, stop_iteration=24, verbose=false)
    run!(simulation)

    @test child.clock.time == parent.clock.time == 2hours
    @test child.clock.iteration == 3 * parent.clock.iteration

    for model in (parent, child)
        @test maximum(abs, model.velocities.u) < 1e-12
        @test maximum(abs, model.velocities.v) < 1e-12
        @test maximum(abs, model.free_surface.displacement) < 1e-12
    end

    T = Array(interior(child.tracers.T))
    @test maximum(abs, T .- T[1:1, 1:1, :]) < 1e-12
end

@testset "A gravity wave crosses the nest" begin
    parent, child, nest = nested_levels(refinement=3)
    bump(λ, φ, z) = 0.1 * exp(-((λ - 2)^2 + (φ - 45)^2) / 2)
    set!(parent, T=stratification, S=35, η=bump)
    set!(child, T=stratification, S=35, η=bump)

    i, j, _ = parent_indices(parent.grid, child.grid)
    block = (nest.refinement, length(i), nest.refinement, length(j))
    amplitude, difference = 0.0, 0.0

    for _ in 1:80
        time_step!(nest, 5minutes)
        η = Array(interior(parent.free_surface.displacement))[i, j, 1]
        blocks = reshape(Array(interior(child.free_surface.displacement)), block)
        restricted_η = sum(blocks, dims=(1, 3))[1, :, 1, :] ./ nest.refinement^2
        amplitude = max(amplitude, maximum(abs, η))
        difference = max(difference, maximum(abs, restricted_η .- η))
    end

    @test difference < 0.05 * amplitude

    # the velocity through the open sides carries the barotropic transport
    u, v = child.velocities
    U, V = child.free_surface.barotropic_velocities
    transport_error = Array(interior(Field(Integral(u, dims=3)))) .- Array(interior(U))
    @test maximum(abs, Array(interior(U))[[1, end], :, :]) > 0.1
    @test maximum(abs, transport_error[[1, end], :, :]) < 1e-12
    transport_error = Array(interior(Field(Integral(v, dims=3)))) .- Array(interior(V))
    @test maximum(abs, transport_error[:, [1, end], :]) < 1e-12
end
