module OnlineNesting

export nested_grid, parent_indices,
       StateExchanger, nested_boundary_conditions, nested_forcings,
       TwoWayNestedModel

using Oceananigans: AbstractModel, prognostic_fields
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: FieldBoundaryConditions, GravityWaveRadiationBoundaryCondition,
                                       LeftBoundary, NormalFlowBoundaryCondition, RightBoundary,
                                       TracerReservoir, ValueBoundaryCondition, regularize_boundary_condition
using Oceananigans.Fields: Field, instantiated_location
using Oceananigans.Grids: Bounded, Center, Face, LatitudeLongitudeGrid, MutableVerticalDiscretization, Open,
                          cpu_face_constructor_r, cpu_face_constructor_x, cpu_face_constructor_y,
                          halo_size, reconstruction_faces
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.TimeSteppers: TimeSteppers, time_step!, update_state!

using NumericalEarth.NestedModels: Interpolated, parent_boundary_conditions, parent_forcings

using DocStringExtensions: TYPEDSIGNATURES

include("nested_grid.jl")
include("state_exchanger.jl")
include("two_way_nested_model.jl")

end # module OnlineNesting
