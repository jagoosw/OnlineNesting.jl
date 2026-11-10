module OnlineNesting

export nested_grid, parent_indices
export EMODnet2024, blend_parent_bottom!, restrict_bottom!

using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: Field, interior, set!
using Oceananigans.Grids: Bounded, LatitudeLongitudeGrid, MutableVerticalDiscretization, Open,
                          cpu_face_constructor_r, cpu_face_constructor_x, cpu_face_constructor_y,
                          Center, halo_size, reconstruction_faces, x_domain, y_domain

using NumericalEarth.DataWrangling: DataWrangling, AbstractStaticBathymetry, Metadatum, bounding_box_suffix,
                                    download_cache, download_with_retries, metadata_path, native_grid
using NumericalEarth.NestedModels: blend_parent_terrain!

using DocStringExtensions: TYPEDSIGNATURES
using Downloads: Downloads
using NCDatasets: NCDataset

include("nested_grid.jl")
include("nested_bathymetry.jl")
include("emodnet.jl")

end # module OnlineNesting
