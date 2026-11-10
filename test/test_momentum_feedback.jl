using Test
using Oceananigans
using Oceananigans.AbstractOperations: Ax, Ay
using Oceananigans.Grids: λnodes, φnodes
using Oceananigans.TimeSteppers: update_state!
using Oceananigans.Units
using OnlineNesting

function nested_levels(parent_grid, grid; width)
    ocean(grid; kwargs...) = HydrostaticFreeSurfaceModel(grid; free_surface = SplitExplicitFreeSurface(grid; substeps=20),
                                                         timestepper = :SplitRungeKutta3,
                                                         tracers = (:T, :S),
                                                         buoyancy = SeawaterBuoyancy(),
                                                         coriolis = HydrostaticSphericalCoriolis(),
                                                         kwargs...)

    parent = ocean(parent_grid)
    exchanger = StateExchanger(parent)
    child = ocean(grid; boundary_conditions = nested_boundary_conditions(exchanger, grid),
                        forcing = nested_forcings(exchanger, grid; rate=1/1hour, width=4))

    return parent, child, TwoWayNestedModel(parent, child, exchanger; width)
end

stratification(λ, φ, z) = 10 + z / 20
surface(model) = Array(interior(model.free_surface.displacement))[:, :, 1]

@testset "The parent takes the nest's transports and wet-mean sea level, $bottom_type" for bottom_type in (GridFittedBottom, PartialCellBottom)
    refinement, width = 3, 2
    parent_grid = LatitudeLongitudeGrid(size=(20, 20, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0))
    underlying_grid = nested_grid(parent_grid; longitude=(3, 7), latitude=(43, 47), depth=50, refinement)

    # land in parts of some parent cells and the whole of others, and a step in the bottom
    bottom(λ, φ) = ifelse(λ > 5.2, 0, ifelse(φ > 45.7, -35, -50))
    grid = ImmersedBoundaryGrid(underlying_grid, bottom_type(bottom))
    parent, child, nest = nested_levels(parent_grid, grid; width)

    wiggle(λ, φ, z) = sind(40λ) * cosd(50φ) * (1 + z / 100)
    set!(parent, T=stratification, S=35, η=-1)
    set!(child, T=stratification, S=35, u=wiggle, v=wiggle, η=wiggle)
    update_state!(nest)

    i, j, k = parent_indices(parent_grid, grid)
    Nx, Ny, Nz = length(i), length(j), length(k)
    inner_i, inner_j = 1+width:Nx-width, 1+width:Ny-width

    x_transport(model) = Array(interior(Field(model.velocities.u * Ax)))
    y_transport(model) = Array(interior(Field(model.velocities.v * Ay)))

    nest_transport = sum(reshape(x_transport(child)[1:refinement:end, :, :], Nx + 1, refinement, Ny, Nz), dims=2)[:, 1, :, :]
    parent_transport = x_transport(parent)[first(i):last(i)+1, j, k]
    @test parent_transport[1+width:Nx-width+1, inner_j, :] ≈ nest_transport[1+width:Nx-width+1, inner_j, :]
    @test all(parent_transport[[1:width; Nx-width+2:Nx+1], :, :] .== 0)

    nest_transport = sum(reshape(y_transport(child)[:, 1:refinement:end, :], refinement, Nx, Ny + 1, Nz), dims=1)[1, :, :, :]
    parent_transport = y_transport(parent)[i, first(j):last(j)+1, k]
    @test parent_transport[inner_i, 1+width:Ny-width+1, :] ≈ nest_transport[inner_i, 1+width:Ny-width+1, :]
    @test all(parent_transport[:, [1:width; Ny-width+2:Ny+1], :] .== 0)

    # the barotropic transport follows the velocity it was given
    U = parent.free_surface.barotropic_velocities.U
    @test Array(interior(U)) ≈ Array(interior(Field(Integral(parent.velocities.u, dims=3))))

    λ, φ = λnodes(grid, Center()), φnodes(grid, Center())
    wet_area = (bottom.(λ, φ') .< 0) .* diff(sind.(φnodes(grid, Face())))'
    block_sum(a) = sum(reshape(a, refinement, Nx, refinement, Ny), dims=(1, 3))[1, :, 1, :]
    wet_mean = block_sum(surface(child) .* wet_area) ./ block_sum(wet_area)
    wet = block_sum(wet_area)[inner_i, inner_j] .> 0

    η = surface(parent)[i, j]
    @test !all(wet) && any(wet)
    @test η[inner_i, inner_j][wet] ≈ wet_mean[inner_i, inner_j][wet]
    @test all(η[inner_i, inner_j][.!wet] .== -1)
end

@testset "A balanced eddy on the nest appears on the parent and stays balanced" begin
    parent_grid = LatitudeLongitudeGrid(size=(40, 40, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0))
    grid = nested_grid(parent_grid; longitude=(3, 7), latitude=(43, 47), refinement=3)
    parent, child, nest = nested_levels(parent_grid, grid; width=2)

    # a geostrophic eddy one parent cell in radius, a tenth of the barotropic deformation radius
    radius = 0.25
    g = parent.free_surface.gravitational_acceleration
    fR(φ) = 2 * parent.coriolis.rotation_rate * sind(φ) * deg2rad(parent_grid.radius)
    η(λ, φ, z=0) = 0.02 * exp(-(((λ - 5) * cosd(45))^2 + (φ - 45)^2) / 2radius^2)
    u(λ, φ, z) =  g / fR(φ) * η(λ, φ) * (φ - 45) / radius^2
    v(λ, φ, z) = -g / fR(φ) * η(λ, φ) * (λ - 5) * cosd(45)^2 / (cosd(φ) * radius^2)

    set!(parent, T=stratification, S=35)
    set!(child; T=stratification, S=35, η, u, v)
    update_state!(nest)

    i, j, _ = parent_indices(parent_grid, grid)
    outside = trues(size(parent_grid, 1), size(parent_grid, 2))
    outside[i, j] .= false

    amplitude = maximum(surface(parent))
    radiated = 0.0

    for _ in 1:48
        time_step!(nest, 5minutes)
        radiated = max(radiated, maximum(abs, surface(parent)[outside]))
    end

    @test amplitude > 0.015
    @test abs(maximum(surface(parent)) - amplitude) < 0.1 * amplitude
    @test radiated < 0.05 * amplitude
end
