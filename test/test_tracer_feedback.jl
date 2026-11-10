using Test
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: λnodes, φnodes, znodes
using Oceananigans.TimeSteppers: update_state!
using OnlineNesting

function ocean(grid; kwargs...)
    return HydrostaticFreeSurfaceModel(grid; free_surface = SplitExplicitFreeSurface(grid; substeps=20),
                                       timestepper = :SplitRungeKutta3,
                                       tracers = (:T, :S, :c),
                                       buoyancy = SeawaterBuoyancy(),
                                       kwargs...)
end

function nested_levels(parent_grid, grid; width)
    parent = ocean(parent_grid)
    exchanger = StateExchanger(parent)
    child = ocean(grid; boundary_conditions = nested_boundary_conditions(exchanger, grid),
                        forcing = nested_forcings(exchanger, grid; rate=1/1hour, width=4))

    return parent, child, TwoWayNestedModel(parent, child, exchanger; width)
end

# Wet cell volumes from the faces of the grid, up to a constant factor
function wet_volumes(grid, wet)
    λ = λnodes(grid, Face(), with_halos=false)
    φ = φnodes(grid, Face(), with_halos=false)
    z = znodes(grid, Face(), with_halos=false)
    Δλ = diff([λ; 2λ[end] - λ[end-1]])[1:size(grid, 1)]
    return wet .* Δλ .* diff(sind.(φ))' .* reshape(diff(z), 1, 1, :)
end

# The volume-weighted mean of the nest's tracer over each parent cell the nest covers
function restricted(tracer, volume, refinement)
    Nx, Ny, Nz = size(volume)
    block = (refinement, Nx ÷ refinement, refinement, Ny ÷ refinement, Nz)
    block_sum(a) = dropdims(sum(reshape(a, block), dims=(1, 3)), dims=(1, 3))
    return block_sum(tracer .* volume) ./ block_sum(volume)
end

stratification(λ, φ, z) = 10 + z / 20

@testset "The parent takes the volume-weighted mean of the nest's wet cells" begin
    refinement, width = 3, 2
    parent_grid = LatitudeLongitudeGrid(size=(20, 20, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0))
    underlying_grid = nested_grid(parent_grid; longitude=(3, 7), latitude=(43, 47), depth=50, refinement)

    # land in parts of some parent cells, the whole of others, and at the bottom of others
    bottom(λ, φ) = ifelse(λ > 5.2, 0, ifelse(φ > 45.7, -30, -50))
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom))
    parent, child, nest = nested_levels(parent_grid, grid; width)

    set!(parent, T=stratification, S=35, c=-1)
    set!(child, T=stratification, S=35, c=(λ, φ, z) -> sind(40λ) * cosd(50φ) * z)
    update_state!(nest)

    λ, φ, z = λnodes(grid, Center()), φnodes(grid, Center()), znodes(grid, Center())
    wet = reshape(z, 1, 1, :) .> bottom.(λ, φ')
    expected = restricted(Array(interior(child.tracers.c)), wet_volumes(grid, wet), refinement)
    wet_fraction = restricted(wet, wet_volumes(grid, trues(size(wet))), refinement)
    wet_parent = wet_fraction .> 0

    i, j, k = parent_indices(parent_grid, grid)
    inner = (1+width:length(i)-width, 1+width:length(j)-width, :)
    @test any(wet_fraction[inner...] .== 0) && any(0 .< wet_fraction[inner...] .< 1) && any(wet_fraction[inner...] .== 1)
    parent_c = Array(interior(parent.tracers.c))
    @test parent_c[i, j, k][inner...][wet_parent[inner...]] ≈ expected[inner...][wet_parent[inner...]]
    @test all(parent_c[i, j, k][inner...][.!wet_parent[inner...]] .== -1)

    parent_c[i[inner[1]], j[inner[2]], k] .= -1
    @test all(parent_c .== -1)
end

@testset "A blob leaves the nest and returns" begin
    refinement, width = 3, 2
    parent_grid = LatitudeLongitudeGrid(size=(24, 16, 4), longitude=(0, 6), latitude=(43, 47), z=(-100, 0),
                                        topology=(Periodic, Bounded, Bounded))
    grid = nested_grid(parent_grid; longitude=(2, 4), latitude=(44, 46), refinement)
    parent, child, nest = nested_levels(parent_grid, grid; width)

    blob(λ, φ, z) = exp(-((λ - 3)^2 + (φ - 45)^2) / (2 * 0.3^2))
    flow(λ, φ, z) = 4 * cosd(φ) # solid-body rotation at 4 m s⁻¹ on the equator
    foreach(model -> set!(model, T=stratification, S=35, c=blob, u=flow), (parent, child))

    volume = wet_volumes(grid, trues(size(grid)))
    i, j, k = parent_indices(parent_grid, grid)
    inner = (1+width:length(i)-width, 1+width:length(j)-width, :)

    Δt = 5minutes
    circuit = round(Int, deg2rad(6) * parent_grid.radius / 4 / Δt)
    restriction_error, least_in_nest = 0.0, 1.0

    for _ in 1:circuit
        time_step!(nest, Δt)
        parent_c, c = Array(interior(parent.tracers.c)), Array(interior(child.tracers.c))
        restriction_error = max(restriction_error, maximum(abs, parent_c[i, j, k][inner...] .- restricted(c, volume, refinement)[inner...]))
        least_in_nest = min(least_in_nest, maximum(c))
    end

    @test restriction_error < 1e-12
    @test least_in_nest < 0.05
    @test maximum(interior(child.tracers.c)) > 0.3
end
