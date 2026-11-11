module OnlineNesting

export nested_grid, nested_timestep, parent_indices,
       StateExchanger, nested_boundary_conditions, nested_forcings,
       TwoWayNestedModel,
       EMODnet2024, refer_to_mean_sea_level!, blend_parent_bottom!, restrict_bottom!,
       TPXO10Atlas, lowest_astronomical_tide, earth_tidal_constituents, earth_tidal_harmonics

using Dates: Dates, DateTime

using Oceananigans: Oceananigans, AbstractModel, TracerTendencyCallsite, fields, prognostic_fields
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Advection: Centered
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.BoundaryConditions: FieldBoundaryConditions, GravityWaveRadiationBoundaryCondition,
                                       LeftBoundary, NormalFlowBoundaryCondition, RightBoundary,
                                       TidalHarmonics, TracerReservoir, ValueBoundaryCondition,
                                       fill_halo_regions!, regularize_boundary_condition
using Oceananigans.Fields: CenterField, Field, OneField, compute!, instantiated_location, interior, interpolate, location, set!
using Oceananigans.Grids: Bounded, Center, Face, Flat, LatitudeLongitudeGrid, MutableVerticalDiscretization, Open, RectilinearGrid,
                          cpu_face_constructor_r, cpu_face_constructor_x, cpu_face_constructor_y,
                          halo_size, inactive_cell, inactive_node, reconstruction_faces, x_domain, y_domain, λnodes, φnodes
using Oceananigans.Models.HydrostaticFreeSurfaceModels: tracer_advecting_velocities
using Oceananigans.MultiRegion: tracer_flux_x, tracer_flux_y
using Oceananigans.Operators: Axᶠᶜᶜ, Ayᶜᶠᶜ, Azᶜᶜᶠ, Vᶜᶜᶜ
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.Simulations: Callback, add_callback!
using Oceananigans.TimeSteppers: TimeSteppers, reconcile_state!, time_step!, update_state!
using Oceananigans.Utils: KernelParameters, launch!

using NumericalEarth.DataWrangling: DataWrangling, AbstractStaticBathymetry, Metadatum, bounding_box_suffix,
                                    download_cache, download_with_retries, inpaint_mask!, metadata_path,
                                    native_grid, unzip
using NumericalEarth.EarthSystemModels: EarthSystemModel
using NumericalEarth.NestedModels: Interpolated, blend_parent_terrain!, parent_boundary_conditions, parent_forcings

using DocStringExtensions: TYPEDSIGNATURES
using Downloads: Downloads
using KernelAbstractions: @index, @kernel
using NCDatasets: NCDataset, defVar

include("nested_grid.jl")
include("nested_timestep.jl")
include("state_exchanger.jl")
include("two_way_nested_model.jl")
include("nested_bathymetry.jl")
include("emodnet.jl")
include("earth_tides.jl")
include("tpxo.jl")

end # module OnlineNesting
