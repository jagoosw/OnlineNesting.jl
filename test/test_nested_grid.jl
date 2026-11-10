using Test
using Oceananigans
using Oceananigans.Grids: Open, rnodes, topology, λnodes, φnodes
using Oceananigans.Operators: Azᶜᶜᶜ
using OnlineNesting

z_faces = [-200, -120, -60, -25, -10, 0]

regular   = LatitudeLongitudeGrid(size=(16, 12, 5), longitude=(-10, -2), latitude=(49, 55), z=z_faces)
stretched = LatitudeLongitudeGrid(size=(6, 12, 5), longitude=[-10, -8.5, -6, -4.5, -3.5, -2.5, -2], latitude=(49, 55), z=z_faces)
mutable   = LatitudeLongitudeGrid(size=(16, 12, 5), longitude=(-10, -2), latitude=(49, 55), z=MutableVerticalDiscretization(z_faces))
immersed  = ImmersedBoundaryGrid(regular, GridFittedBottom((λ, φ) -> -150))

parent_grids = (; regular, stretched, mutable, immersed)

@testset "nested_grid aligned with a $name parent, refinement $refinement" for name in keys(parent_grids), refinement in (2, 5)
    parent_grid = parent_grids[name]
    grid = nested_grid(parent_grid; longitude=(-7.1, -3.4), latitude=(50.4, 53.1), refinement, depth=70)
    i, j, k = parent_indices(parent_grid, grid)

    @test size(grid) == (refinement * length(i), refinement * length(j), length(k))
    @test topology(grid, 1) <: Oceananigans.Grids.TwoSided{Open, Open}
    @test (grid.z isa MutableVerticalDiscretization) == (name == :mutable)

    faces(nodes, grid) = collect(nodes(grid, Face()))
    @test faces(λnodes, grid)[1:refinement:end] ≈ faces(λnodes, parent_grid)[first(i):last(i)+1]
    @test faces(φnodes, grid)[1:refinement:end] ≈ faces(φnodes, parent_grid)[first(j):last(j)+1]
    @test faces(rnodes, grid) == faces(rnodes, parent_grid)[first(k):last(k)+1]

    # the nest tiles its parent: the cells within each parent cell have its area
    Nx, Ny, _ = size(grid)
    nest_area = [Azᶜᶜᶜ(m, n, 1, grid) for m in 1:Nx, n in 1:Ny]
    tiled_area = [sum(nest_area[refinement*(m-1)+1:refinement*m, refinement*(n-1)+1:refinement*n]) for m in 1:length(i), n in 1:length(j)]
    parent_area = [Azᶜᶜᶜ(m, n, 1, parent_grid) for m in i, n in j]
    @test tiled_area ≈ parent_area
end
