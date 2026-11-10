module OnlineNesting

export nested_grid, parent_indices

using Oceananigans.Architectures: architecture
using Oceananigans.Grids: Bounded, LatitudeLongitudeGrid, MutableVerticalDiscretization, Open,
                          cpu_face_constructor_r, cpu_face_constructor_x, cpu_face_constructor_y,
                          halo_size, reconstruction_faces

using DocStringExtensions: TYPEDSIGNATURES

include("nested_grid.jl")

end # module OnlineNesting
