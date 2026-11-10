module OnlineNesting

export nested_grid, parent_indices

using Oceananigans.Architectures: CPU, architecture, on_architecture
using Oceananigans.Grids: Bounded, Face, LatitudeLongitudeGrid, MutableVerticalDiscretization, Open,
                          halo_size, rnodes, λnodes, φnodes

using DocStringExtensions: TYPEDSIGNATURES

include("nested_grid.jl")

end # module OnlineNesting
