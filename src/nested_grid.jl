nearest_face(faces, x) = argmin(abs.(faces .- x))

# A regularly spaced dimension is constructed from its extent, a stretched one from every face
refined_faces(::Tuple, faces, refinement) = (first(faces), last(faces))

function refined_faces(constructor, faces, refinement)
    fractions = (0:refinement-1) ./ refinement
    refined = [faces[n] + (faces[n+1] - faces[n]) * f for n in 1:length(faces)-1 for f in fractions]
    return push!(refined, last(faces))
end

vertical_coordinate(::MutableVerticalDiscretization, faces) = MutableVerticalDiscretization(faces)
vertical_coordinate(coordinate, faces) = faces

function cpu_faces(grid)
    Nx, Ny, Nz = size(grid)
    return (reconstruction_faces(cpu_face_constructor_x(grid), Nx),
            reconstruction_faces(cpu_face_constructor_y(grid), Ny),
            reconstruction_faces(cpu_face_constructor_r(grid), Nz))
end

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
nested_grid(parent_grid; longitude=(-6.2, -3.9), latitude=(51, 53), refinement=5, depth=100)

# output
10×10×2 LatitudeLongitudeGrid{Float64, Open, Open, Bounded} on CPU with 3×3×3 halo
├── longitude: Open     λ ∈ [-6.0, -4.0]  regularly spaced with Δλ=0.2
├── latitude:  Open     φ ∈ [51.0, 53.0]  regularly spaced with Δφ=0.2
└── z:         Bounded  z ∈ [-100.0, 0.0] regularly spaced with Δz=50.0
```
"""
function nested_grid(parent_grid; longitude, latitude, refinement, depth = parent_grid.Lz)
    λᶠ, φᶠ, zᶠ = cpu_faces(parent_grid)

    i = nearest_face(λᶠ, longitude[1]):nearest_face(λᶠ, longitude[2])
    j = nearest_face(φᶠ, latitude[1]):nearest_face(φᶠ, latitude[2])
    k = nearest_face(zᶠ, -depth):length(zᶠ)

    size = (refinement * (length(i) - 1), refinement * (length(j) - 1), length(k) - 1)
    z = refined_faces(cpu_face_constructor_r(parent_grid), zᶠ[k], 1)

    return LatitudeLongitudeGrid(architecture(parent_grid), eltype(parent_grid);
                                 size,
                                 longitude = refined_faces(cpu_face_constructor_x(parent_grid), λᶠ[i], refinement),
                                 latitude = refined_faces(cpu_face_constructor_y(parent_grid), φᶠ[j], refinement),
                                 z = vertical_coordinate(parent_grid.z, z),
                                 halo = halo_size(parent_grid),
                                 radius = parent_grid.radius,
                                 topology = (Open, Open, Bounded))
end

"""
$(TYPEDSIGNATURES)

Return the ranges `(i, j, k)` of the cells of `parent_grid` that `nest_grid` covers.
"""
function parent_indices(parent_grid, nest_grid)
    return map(cpu_faces(parent_grid), cpu_faces(nest_grid)) do parent_faces, nest_faces
        nearest_face(parent_faces, first(nest_faces)):nearest_face(parent_faces, last(nest_faces))-1
    end
end
