using Test
using Oceananigans
using Oceananigans.Grids: Open, TwoSided, halo_size, rnodes, topology, λnodes, φnodes
using OnlineNesting

z_faces = [-200, -120, -60, -25, -10, 0]
stretched_longitude = [-10, -8.5, -6, -4.5, -3.5, -2.5, -2]

regular   = LatitudeLongitudeGrid(size=(16, 12, 5), longitude=(-10, -2), latitude=(49, 55), z=z_faces)
stretched = LatitudeLongitudeGrid(size=(6, 12, 5), longitude=stretched_longitude, latitude=(49, 55), z=z_faces, halo=(4, 4, 4))
mutable   = LatitudeLongitudeGrid(size=(16, 12, 5), longitude=(-10, -2), latitude=(49, 55), z=MutableVerticalDiscretization(z_faces))
immersed  = ImmersedBoundaryGrid(regular, GridFittedBottom((λ, φ) -> -150))
single    = LatitudeLongitudeGrid(Float32; size=(16, 12, 5), longitude=(-10, -2), latitude=(49, 55), z=z_faces, radius=1f6)

parent_grids = (; regular, stretched, mutable, immersed, single)

faces(nodes, grid) = collect(nodes(grid, Face()))
nearest(values, x) = values[argmin(abs.(values .- x))]

@testset "nested_grid in a $name parent, refinement $refinement" for name in keys(parent_grids), refinement in (2, 5)
    parent_grid = parent_grids[name]
    longitude, latitude, depth = (-7.1, -3.4), (50.4, 53.1), 70
    grid = nested_grid(parent_grid; longitude, latitude, refinement, depth)
    i, j, k = parent_indices(parent_grid, grid)

    # the nest spans the parent faces nearest to the requested extents
    for (nodes, extent) in ((λnodes, longitude), (φnodes, latitude), (rnodes, (-depth, 0)))
        parent_faces, nest_faces = faces(nodes, parent_grid), faces(nodes, grid)
        @test first(nest_faces) ≈ nearest(parent_faces, extent[1])
        @test last(nest_faces) ≈ nearest(parent_faces, extent[2])
    end

    # and splits each parent cell it covers into `refinement` equal cells
    @test diff(faces(λnodes, grid)) ≈ repeat(diff(faces(λnodes, parent_grid))[i] ./ refinement, inner=refinement)
    @test diff(faces(φnodes, grid)) ≈ repeat(diff(faces(φnodes, parent_grid))[j] ./ refinement, inner=refinement)
    @test diff(faces(rnodes, grid)) == diff(faces(rnodes, parent_grid))[k]

    @test topology(grid) == (TwoSided{Open, Open}, TwoSided{Open, Open}, Bounded)
    @test (grid.z isa MutableVerticalDiscretization) == (name == :mutable)
    @test (grid.Δλᶠᵃᵃ isa Number) == (name != :stretched)
    @test (halo_size(grid), eltype(grid), grid.radius) == (halo_size(parent_grid), eltype(parent_grid), parent_grid.radius)
end
