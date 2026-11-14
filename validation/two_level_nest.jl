# Two-level nested test: shelf → Bristol Channel on GPU, no Severn, no land model.
# If this crashes, the basic nesting machinery is the problem.
# If it passes, the three-level depth or Severn-specific setup is the trigger.
#
#   julia --project=test -t auto validation/two_level_nest.jl

using Dates
using CDSAPI
using CopernicusMarine
using Oceananigans
using Oceananigans.Advection: WENO
using Oceananigans.Units
using Oceananigans.BoundaryConditions: TracerReservoir, tidal_boundary_conditions
using Oceananigans.Fields: interpolate
using Oceananigans.Grids: Open, MutableVerticalDiscretization
using NumericalEarth
using NumericalEarth.DataWrangling: BoundingBox, GLORYSDaily, Metadata, MetadataSet
using NumericalEarth.NestedModels: Interpolated
using OnlineNesting
using Printf

using NumericalEarth.EarthSystemModels.InterfaceComputations: ReynoldsScalingFunction,
    MomentumRoughnessLength, gravity_wave_parameter, compute_air_kinematic_viscosity
import NumericalEarth.EarthSystemModels.InterfaceComputations: roughness_length

@inline function roughness_length(ℓ::MomentumRoughnessLength{FT}, u★, Uᵃᵗ, ℂᵃᵗ=nothing, Tᵃᵗ=nothing) where FT
    u★_pos = max(u★, zero(u★))
    ν = compute_air_kinematic_viscosity(ℓ.air_kinematic_viscosity, ℂᵃᵗ, Tᵃᵗ)
    g = ℓ.gravitational_acceleration
    ℂg = gravity_wave_parameter(ℓ.wave_formulation, Uᵃᵗ)
    ℂν = ℓ.smooth_wall_parameter
    ℓᵂ = ℂg * u★_pos^2 / g
    ℓᴿ = ifelse(u★_pos == zero(u★), zero(u★), ℂν * ν / u★_pos)
    ℓ★ = ℓᵂ + ℓᴿ
    ℓ_max = ℓ.maximum_roughness_length
    return min(ℓ★, ℓ_max)
end

@inline (s::ReynoldsScalingFunction)(R★, args...) =
    ifelse(R★ <= 0, convert(eltype(R★), 0), s.A / max(R★, eps(eltype(R★))) ^ s.b)

#####
##### Configuration
#####

arch = GPU()
start_date = DateTime(2018, 1, 1)
end_date   = DateTime(2021, 1, 1)

tpxo_dir    = get(ENV, "TPXO_DIR", "/home/js5256/projects/OnlineNesting.jl/data/TPXO10_atlas_v2_nc/")
glorys_dir  = get(ENV, "GLORYS_DIR", "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/GLORYS/")
era5_dir    = get(ENV, "ERA5_DIR", "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/ERA5/")
emodnet_dir = get(ENV, "EMODNET_DIR", "/home/js5256/projects/OnlineNesting.jl/data/EMODnet/")

refinement    = 5
width         = 3
minimum_depth = 10
Nz            = 50

#####
##### Vertical grid
#####

z_faces = [-250 * (1 - k / Nz)^1.8 for k in 0:Nz]
z = MutableVerticalDiscretization(z_faces)

#####
##### Grids (shelf + Bristol only)
#####

shelf_grid = LatitudeLongitudeGrid(arch;
    size = (340, 320, Nz),
    longitude = (-12.0, 5.0),
    latitude = (46.0, 62.0),
    halo = (7, 7, 7),
    topology = (Open, Open, Bounded),
    z)

bristol_grid = nested_grid(shelf_grid;
    longitude = (-5.5, -2.0),
    latitude = (50.8, 52.2),
    depth = 250,
    refinement)

#####
##### Timestep
#####

result = nested_timestep((shelf_grid, bristol_grid), (0.5, 1.0))
Δt = result.timesteps[1]

#####
##### Bathymetry
#####

tpxo = TPXO10Atlas(; dir=tpxo_dir)
lat_field = lowest_astronomical_tide(tpxo, shelf_grid)
λ_west = Oceananigans.Grids.x_domain(lat_field.grid)[1]
lowest_tide(λ, φ) = interpolate((mod(λ - λ_west, 360) + λ_west, φ), lat_field)

function emodnet_bottom_height(grid; interpolation_passes, major_basins=Inf, lat_correction=lowest_tide, dir=emodnet_dir)
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(grid), dir)
    bottom_height = regrid_bathymetry(grid, metadatum; interpolation_passes, major_basins)
    return refer_to_mean_sea_level!(bottom_height, lat_correction; minimum_depth)
end

shelf_bottom   = emodnet_bottom_height(shelf_grid;   interpolation_passes=4, major_basins=1, lat_correction=0)
bristol_bottom = emodnet_bottom_height(bristol_grid;  interpolation_passes=2)

restrict_bottom!(shelf_bottom, bristol_bottom; width, minimum_depth)

immersed(grid, bottom) = ImmersedBoundaryGrid(grid, PartialCellBottom(Array(interior(bottom, :, :, 1))))
shelf_grid, bristol_grid = immersed.((shelf_grid, bristol_grid), (shelf_bottom, bristol_bottom))

#####
##### Tidal forcing on the shelf boundary
#####

harmonics  = earth_tidal_harmonics(start_date; ramp_time=6hours)
tidal_bcs  = tidal_boundary_conditions(shelf_grid, harmonics, tpxo)

#####
##### GLORYS open boundary conditions on the shelf
#####

glorys = GLORYSDaily()
glorys_region = BoundingBox(longitude = (-12.5, 5.5), latitude = (45.5, 62.5))
glorys_kw = (; dataset=glorys, start_date, end_date, region=glorys_region, dir=glorys_dir)
fts_kw = (; time_indices_in_memory=30)

glorys_T = FieldTimeSeries(Metadata(:temperature; glorys_kw...), arch; fts_kw...)
glorys_S = FieldTimeSeries(Metadata(:salinity;    glorys_kw...), arch; fts_kw...)
glorys_u = FieldTimeSeries(Metadata(:u_velocity;  glorys_kw...), arch; fts_kw...)
glorys_v = FieldTimeSeries(Metadata(:v_velocity;  glorys_kw...), arch; fts_kw...)

parent(glorys_S.data) .= max.(parent(glorys_S.data), 0)

reservoir(source) = ValueBoundaryCondition(Interpolated(source); scheme=TracerReservoir())
normal(source)    = NormalFlowBoundaryCondition(Interpolated(source))

shelf_boundary_conditions = (
    U = tidal_bcs.U,
    V = tidal_bcs.V,
    η = tidal_bcs.η,
    u = FieldBoundaryConditions(west=normal(glorys_u), east=normal(glorys_u)),
    v = FieldBoundaryConditions(south=normal(glorys_v), north=normal(glorys_v)),
    T = FieldBoundaryConditions(west=reservoir(glorys_T), east=reservoir(glorys_T),
                                south=reservoir(glorys_T), north=reservoir(glorys_T)),
    S = FieldBoundaryConditions(west=reservoir(glorys_S), east=reservoir(glorys_S),
                                south=reservoir(glorys_S), north=reservoir(glorys_S)),
)

#####
##### Ocean simulations
#####

tracers = (:T, :S)
vi = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)
tracer_advection = WENO(order=7, time_discretization=vi, minimum_buffer_upwind_order=1)

shelf_ocean = ocean_simulation(shelf_grid; tracers, tracer_advection,
    boundary_conditions=shelf_boundary_conditions, warn=false)

glorys_init_kw = (; dataset=glorys, date=start_date, region=glorys_region, dir=glorys_dir)
set!(shelf_ocean.model, T = Metadatum(:temperature; glorys_init_kw...),
                        S = Metadatum(:salinity;    glorys_init_kw...))

#####
##### ERA5 atmosphere and radiation
#####

era5 = ERA5YearlySingleLevel()
era5_kw = (; dataset=era5, dir=era5_dir)

shelf_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date, pressure=:mean_sea_level_pressure, era5_kw...)
shelf_radiation  = ERA5PrescribedRadiation(arch;  start_date, end_date, era5_kw...)

bristol_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date, pressure=:mean_sea_level_pressure, era5_kw...)
bristol_radiation  = ERA5PrescribedRadiation(arch;  start_date, end_date, era5_kw...)

#####
##### EarthSystemModels
#####

shelf_model = OceanOnlyModel(shelf_ocean;
    atmosphere=shelf_atmosphere, radiation=shelf_radiation)
parent(shelf_ocean.model.tracers.S) .= max.(parent(shelf_ocean.model.tracers.S), 0)

bristol_exchanger = StateExchanger(shelf_model)
bristol_sponge    = nested_forcings(bristol_exchanger, bristol_grid;
    rate=1/1hours, width=refinement*width)
bristol_bcs       = nested_boundary_conditions(bristol_exchanger, bristol_grid)

bristol_ocean = ocean_simulation(bristol_grid; tracers, tracer_advection,
    boundary_conditions=bristol_bcs, forcing=bristol_sponge, warn=false)
set!(bristol_ocean.model, T = Metadatum(:temperature; glorys_init_kw...),
                          S = Metadatum(:salinity;    glorys_init_kw...))

bristol_model = OceanOnlyModel(bristol_ocean;
    atmosphere=bristol_atmosphere, radiation=bristol_radiation)
parent(bristol_ocean.model.tracers.S) .= max.(parent(bristol_ocean.model.tracers.S), 0)

#####
##### Two-level nested model
#####

nest = TwoWayNestedModel(shelf_model, bristol_model, bristol_exchanger; width)

#####
##### One time step
#####

@printf("Two-level nest test: Δt = %.1f s\n", Δt)
@printf("  shelf:   %3d × %3d × %d\n", size(shelf_grid)...)
@printf("  bristol: %3d × %3d × %d\n", size(bristol_grid)...)

@info "Starting single time step on GPU..."
time_step!(nest, Δt)
@info "Time step completed without error."
