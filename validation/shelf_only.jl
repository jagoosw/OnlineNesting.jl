# Shelf-only test: run the outermost shelf model with TPXO tidal boundaries for one
# time step on GPU, without any nesting. If this crashes, the bug is in the boundary
# conditions or ocean setup, not in OnlineNesting's two-way nesting.
#
#   julia --project=test -t auto validation/shelf_only.jl

using Dates
using CDSAPI
using CopernicusMarine
using Oceananigans
using Oceananigans.Advection: WENO
using Oceananigans.Units
using Oceananigans.BoundaryConditions: TracerReservoir, tidal_boundary_conditions
using Oceananigans.Grids: Open, MutableVerticalDiscretization
using NumericalEarth
using NumericalEarth.DataWrangling: BoundingBox, GLORYSDaily, Metadata, MetadataSet
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
dates = (start_date, end_date)

tpxo_dir    = get(ENV, "TPXO_DIR", "/home/js5256/projects/OnlineNesting.jl/data/TPXO10_atlas_v2_nc/")
glorys_dir  = get(ENV, "GLORYS_DIR", "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/GLORYS/")
era5_dir    = get(ENV, "ERA5_DIR", "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/ERA5/")
emodnet_dir = get(ENV, "EMODNET_DIR", "/home/js5256/projects/OnlineNesting.jl/data/EMODnet/")

minimum_depth = 10
Nz = 50

#####
##### Vertical grid: 50 surface-intensified layers to 250 m
#####

z_faces = [-250 * (1 - k / Nz)^1.8 for k in 0:Nz]
z = MutableVerticalDiscretization(z_faces)

#####
##### Shelf grid
#####

shelf_grid = LatitudeLongitudeGrid(arch;
    size = (340, 320, Nz),
    longitude = (-12.0, 5.0),
    latitude = (46.0, 62.0),
    halo = (7, 7, 7),
    topology = (Open, Open, Bounded),
    z)

#####
##### Bathymetry (EMODnet referenced to MSL via TPXO lowest astronomical tide)
#####

tpxo = TPXO10Atlas(; dir=tpxo_dir)

using Oceananigans.Fields: interpolate
lat_field = lowest_astronomical_tide(tpxo, shelf_grid)
λ_west = Oceananigans.Grids.x_domain(lat_field.grid)[1]
lowest_tide(λ, φ) = interpolate((mod(λ - λ_west, 360) + λ_west, φ), lat_field)

metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(shelf_grid), dir=emodnet_dir)
shelf_bottom = regrid_bathymetry(shelf_grid, metadatum; interpolation_passes=4, major_basins=1)
refer_to_mean_sea_level!(shelf_bottom, 0; minimum_depth)

shelf_grid = ImmersedBoundaryGrid(shelf_grid, PartialCellBottom(Array(interior(shelf_bottom, :, :, 1))))

#####
##### Tidal forcing on the shelf boundary
#####

harmonics = earth_tidal_harmonics(start_date; ramp_time=6hours)
tidal_bcs = tidal_boundary_conditions(shelf_grid, harmonics, tpxo)

#####
##### GLORYS reanalysis open boundary conditions on the shelf
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

using NumericalEarth.NestedModels: Interpolated
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
##### Ocean simulation (CATKE closure, WENO advection, SplitRungeKutta3)
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

#####
##### EarthSystemModel (OceanOnlyModel wraps the ocean with atmosphere + radiation)
#####

shelf_model = OceanOnlyModel(shelf_ocean;
    atmosphere=shelf_atmosphere, radiation=shelf_radiation)

parent(shelf_ocean.model.tracers.S) .= max.(parent(shelf_ocean.model.tracers.S), 0)

#####
##### Timestep from CFL
#####

Δx_min = Oceananigans.Grids.minimum_xspacing(shelf_grid)
Δy_min = Oceananigans.Grids.minimum_yspacing(shelf_grid)
Δ_min  = min(Δx_min, Δy_min)
c_max  = sqrt(9.81 * 250)
Δt     = 0.5 * Δ_min / c_max
@printf("Shelf-only test: Δt = %.1f s, grid = %d × %d × %d\n", Δt, size(shelf_grid)...)

#####
##### One time step
#####

@info "Starting single time step on GPU..."
time_step!(shelf_model, Δt)
@info "Time step completed without error."
