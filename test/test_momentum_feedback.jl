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

@testset "The parent takes the nest's velocity and wet-mean sea level, $bottom_type" for bottom_type in (GridFittedBottom, PartialCellBottom)
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
    inner = (1+width:Nx-width, 1+width:Ny-width)

    # the nest's transport through each parent face, over the area of the nest's faces there
    function restricted(velocity, area, on_face, dims)
        transport = Array(interior(Field(velocity * area)))[on_face..., :]
        face_area = Array(interior(Field(velocity * 0 + area)))[on_face..., :]
        block = dims == 2 ? (Nx + 1, refinement, Ny, Nz) : (refinement, Nx, Ny + 1, Nz)
        return dropdims(sum(reshape(transport, block); dims) ./ sum(reshape(face_area, block); dims); dims)
    end

    faces = 1+width:Nx-width+1
    nest_u = restricted(child.velocities.u, Ax, (1:refinement:size(grid, 1)+1, :), 2)
    parent_u = Array(interior(parent.velocities.u))[first(i):last(i)+1, j, k]
    @test maximum(abs, nest_u) > 0.1
    @test parent_u[faces, inner[2], :] ≈ nest_u[faces, inner[2], :]
    @test all(parent_u[[1:width; Nx-width+2:Nx+1], :, :] .== 0)

    faces = 1+width:Ny-width+1
    nest_v = restricted(child.velocities.v, Ay, (:, 1:refinement:size(grid, 2)+1), 1)
    parent_v = Array(interior(parent.velocities.v))[i, first(j):last(j)+1, k]
    @test parent_v[inner[1], faces, :] ≈ nest_v[inner[1], faces, :]
    @test all(parent_v[:, [1:width; Ny-width+2:Ny+1], :] .== 0)

    # the barotropic transport follows the velocity it was given
    U = parent.free_surface.barotropic_velocities.U
    @test Array(interior(U)) ≈ Array(interior(Field(Integral(parent.velocities.u, dims=3))))

    λ, φ = λnodes(grid, Center()), φnodes(grid, Center())
    wet_area = (bottom.(λ, φ') .< 0) .* diff(sind.(φnodes(grid, Face())))'
    block_sum(a) = sum(reshape(a, refinement, Nx, refinement, Ny), dims=(1, 3))[1, :, 1, :]
    wet_mean = block_sum(surface(child) .* wet_area) ./ block_sum(wet_area)
    wet = block_sum(wet_area)[inner...] .> 0

    η = surface(parent)[i, j]
    @test !all(wet) && any(wet)
    @test η[inner...][wet] ≈ wet_mean[inner...][wet]
    @test all(η[inner...][.!wet] .== -1)
end

@testset "A parent face thinner than the nest's takes no more than the nest's velocity" begin
    underlying_grid = LatitudeLongitudeGrid(size=(20, 20, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0))
    thin_cells = PartialCellBottom((λ, φ) -> ifelse(λ > 5, -75.01, -100), minimum_fractional_cell_height=0)
    parent_grid = ImmersedBoundaryGrid(underlying_grid, thin_cells)
    grid = nested_grid(underlying_grid; longitude=(3, 7), latitude=(43, 47), refinement=3)
    parent, child, nest = nested_levels(parent_grid, grid; width=2)

    set!(parent, T=stratification, S=35)
    set!(child, T=stratification, S=35, u=0.1, v=0.1)
    update_state!(nest)

    @test 0.099 < maximum(parent.velocities.u) < 0.101
    @test 0.099 < maximum(parent.velocities.v) < 0.101
end

# The largest sea level left on the parent, relative to its start, when it runs alone from a nest's eddy
function remaining_eddy(; balanced)
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
    balanced || set!(parent, u=0, v=0)

    amplitude = maximum(surface(parent))
    remaining = amplitude

    for _ in 1:48
        time_step!(parent, 5minutes)
        remaining = min(remaining, maximum(surface(parent)))
    end

    return remaining / amplitude
end

@testset "A balanced eddy on the nest is balanced on the parent" begin
    @test remaining_eddy(balanced=true) > 0.75
    @test remaining_eddy(balanced=false) < 0.25
end
