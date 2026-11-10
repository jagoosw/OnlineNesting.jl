nearest_face(faces, x) = argmin(abs.(faces .- x))

# Faces `indices` of `faces` with each cell split into `refinement` equal parts: the extent of
# regularly spaced faces, which the grid constructor divides itself, and every face otherwise.
refined_faces(faces::AbstractRange, indices, refinement) = (faces[first(indices)], faces[last(indices)])

function refined_faces(faces, indices, refinement)
    cells = first(indices):last(indices)-1
    fractions = (0:refinement-1) ./ refinement
    refined = [faces[i] + (faces[i+1] - faces[i]) * f for i in cells for f in fractions]
    return push!(refined, faces[last(indices)])
end

vertical_coordinate(::MutableVerticalDiscretization, faces) = MutableVerticalDiscretization(faces)
vertical_coordinate(coordinate, faces) = faces

cpu_faces(nodes, grid) = on_architecture(CPU(), nodes(grid, Face()))

"""
$(TYPEDSIGNATURES)

Return a `LatitudeLongitudeGrid` with open sides nested in `parent_grid`, with `refinement` cells
per parent cell in each horizontal direction. The `longitude` and `latitude` extents and the
`depth` [m] are moved to the nearest parent faces, so every parent face within the nest is
also a nest face. The nest has the parent's vertical faces above `depth`, and its halo, radius
and kind of vertical coordinate.

```jldoctest
using Oceananigans
using OnlineNesting

parent_grid = LatitudeLongitudeGrid(size=(8, 6, 4), longitude=(-10, -2), latitude=(49, 55), z=(-200, 0))
grid = nested_grid(parent_grid; longitude=(-6.2, -3.9), latitude=(51, 53), refinement=5, depth=100)
size(grid)

# output
(10, 10, 2)
```
"""
function nested_grid(parent_grid; longitude, latitude, refinement, depth = parent_grid.Lz)
    λᶠ = cpu_faces(λnodes, parent_grid)
    φᶠ = cpu_faces(φnodes, parent_grid)
    zᶠ = cpu_faces(rnodes, parent_grid)

    i = nearest_face(λᶠ, longitude[1]):nearest_face(λᶠ, longitude[2])
    j = nearest_face(φᶠ, latitude[1]):nearest_face(φᶠ, latitude[2])
    k = nearest_face(zᶠ, -depth):length(zᶠ)

    size = (refinement * (length(i) - 1), refinement * (length(j) - 1), length(k) - 1)
    z = vertical_coordinate(parent_grid.z, refined_faces(zᶠ, k, 1))

    return LatitudeLongitudeGrid(architecture(parent_grid), eltype(parent_grid);
                                 size, z,
                                 longitude = refined_faces(λᶠ, i, refinement),
                                 latitude = refined_faces(φᶠ, j, refinement),
                                 halo = halo_size(parent_grid),
                                 radius = parent_grid.radius,
                                 topology = (Open, Open, Bounded))
end

"""
$(TYPEDSIGNATURES)

Return the ranges `(i, j, k)` of the cells of `parent_grid` that `nest_grid` covers.
"""
function parent_indices(parent_grid, nest_grid)
    return map((λnodes, φnodes, rnodes)) do nodes
        parent_faces = cpu_faces(nodes, parent_grid)
        nest_faces = cpu_faces(nodes, nest_grid)
        nearest_face(parent_faces, first(nest_faces)):nearest_face(parent_faces, last(nest_faces))-1
    end
end
