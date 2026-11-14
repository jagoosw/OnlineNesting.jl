# Three-level nested simulation: NW European shelf → Bristol Channel → Severn Estuary.
#
# The innermost domain resolves the Severn Estuary at ~140 m from Llantwit Major
# to the tidal limit east of Gloucester. The middle domain covers the Bristol
# Channel at ~690 m and the outer domain the NW European shelf at ~3.5 km.
#
# Prerequisites:
#   - TPXO10 atlas files in TPXO_DIR (licensed; request at https://www.tpxo.net)
#   - Copernicus Marine credentials: COPERNICUSMARINE_SERVICE_USERNAME, _PASSWORD
#   - ERA5 data auto-downloads; set ERA5_DIR for a shared cache
#   - GloFAS river routing (optional): set GLOFAS_DIR to a directory with GloFAS reanalysis data
#
#   julia --project=test -t auto validation/severn_estuary.jl [days, default 30]
#
# For GPU, load CUDA before this script and set arch = GPU() below.

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

# TODO: fix upstream in NumericalEarth. The MOST iteration can produce negative u★ when the
# similarity profile flips sign under extreme atmosphere-ocean temperature differences at
# startup. A negative u★ makes the viscous roughness term ℂν ν / u★ large and negative,
# so the momentum roughness length ℓu goes negative and log(h / ℓu) throws DomainError.
# Even when ℓu stays positive, R★ = ℓu u★ / ν is negative and R★^0.72 also throws.
# Guard both by clamping u★ ≥ 0 in the roughness length; the iteration then converges to
# zero fluxes for those cells.
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

arch = CPU()
days_to_run = isempty(ARGS) ? 30 : parse(Float64, ARGS[1])
start_date = DateTime(2018, 1, 1)
end_date   = DateTime(2021, 1, 1)
dates = (start_date, end_date)

output_dir  = "/nfs/roberts/scratch/pi_ey239/js5256/severn"
tpxo_dir    = "/home/js5256/projects/OnlineNesting.jl/data/TPXO10_atlas_v2_nc/"
glorys_dir = "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/GLORYS/"
era5_dir    = "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/ERA5/"
glofas_dir  = "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/GloFAS/"
emodnet_dir = "/home/js5256/projects/OnlineNesting.jl/data/EMODnet/"
mkpath(output_dir)

refinement    = 5
width         = 3
minimum_depth = 10
Nz            = 50

#####
##### Vertical grid: 50 surface-intensified layers to 250 m
#####

z_faces = [-250 * (1 - k / Nz)^1.8 for k in 0:Nz]
z = MutableVerticalDiscretization(z_faces)

#####
##### Horizontal grids
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

severn_grid = nested_grid(bristol_grid;
    longitude = (-3.5, -2.5),
    latitude = (51.1, 51.75),
    depth = 100,
    refinement)

#####
##### Timestep from CFL constraints (before immersing the grids)
#####

result = nested_timestep(
    (shelf_grid, bristol_grid, severn_grid),
    (0.5, 1.0, 2.0))
Δt = result.timesteps[1] / 4

#####
##### Bathymetry (EMODnet referenced to MSL via TPXO lowest astronomical tide)
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
severn_bottom  = emodnet_bottom_height(severn_grid;   interpolation_passes=1)

blend_parent_bottom!(severn_bottom,  bristol_bottom; width, minimum_depth)
restrict_bottom!(bristol_bottom, severn_bottom;  width, minimum_depth)
restrict_bottom!(shelf_bottom,   bristol_bottom; width, minimum_depth)

sea = map(b -> Array(interior(b, :, :, 1)) .< 0,
          (shelf = shelf_bottom, bristol = bristol_bottom, severn = severn_bottom))

immersed(grid, bottom) = ImmersedBoundaryGrid(grid, PartialCellBottom(Array(interior(bottom, :, :, 1))))
shelf_grid, bristol_grid, severn_grid =
    immersed.((shelf_grid, bristol_grid, severn_grid),
              (shelf_bottom, bristol_bottom, severn_bottom))

#####
##### Tidal forcing on the shelf boundary
#####

harmonics  = earth_tidal_harmonics(start_date; ramp_time=6hours)
tidal_bcs  = tidal_boundary_conditions(shelf_grid, harmonics, tpxo)

#####
##### GLORYS reanalysis open boundary conditions on the shelf
#####

glorys = GLORYSDaily()
glorys_region = BoundingBox(longitude = (-12.5, 5.5), latitude = (45.5, 62.5))
glorys_kw = (; dataset=glorys, start_date, end_date, region=glorys_region)
isnothing(glorys_dir) || (glorys_kw = merge(glorys_kw, (; dir=glorys_dir)))
fts_kw = (; time_indices_in_memory=30)

glorys_T = FieldTimeSeries(Metadata(:temperature; glorys_kw...), arch; fts_kw...)
glorys_S = FieldTimeSeries(Metadata(:salinity;    glorys_kw...), arch; fts_kw...)
glorys_u = FieldTimeSeries(Metadata(:u_velocity;  glorys_kw...), arch; fts_kw...)
glorys_v = FieldTimeSeries(Metadata(:v_velocity;  glorys_kw...), arch; fts_kw...)

# GLORYS reanalysis can carry slightly negative salinity near coasts; TEOS10 requires S_A ≥ 0
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
##### ERA5 atmosphere and radiation per level
#####

era5 = ERA5YearlySingleLevel()
era5_kw = isnothing(era5_dir) ? (; dataset=era5) : (; dataset=era5, dir=era5_dir)

shelf_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date, pressure=:mean_sea_level_pressure, era5_kw...)
shelf_radiation  = ERA5PrescribedRadiation(arch;  start_date, end_date, era5_kw...)

bristol_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date, pressure=:mean_sea_level_pressure, era5_kw...)
bristol_radiation  = ERA5PrescribedRadiation(arch;  start_date, end_date, era5_kw...)

severn_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date, pressure=:mean_sea_level_pressure, era5_kw...)
severn_radiation  = ERA5PrescribedRadiation(arch;  start_date, end_date, era5_kw...)

#####
##### GloFAS river routing
#####

shelf_land = GloFASPrescribedLand(shelf_grid;
                                  start_date, end_date,
                                  region=BoundingBox(shelf_grid), dir=glofas_dir)
shelf_river_routing = shelf_land.river_routing

bristol_land = GloFASPrescribedLand(bristol_grid;
                                    start_date, end_date,
                                    region=BoundingBox(bristol_grid), dir=glofas_dir)
bristol_river_routing = bristol_land.river_routing

severn_land = GloFASPrescribedLand(severn_grid;
                                   start_date, end_date,
                                   region=BoundingBox(severn_grid), dir=glofas_dir)
severn_river_routing = severn_land.river_routing

#####
##### Ocean simulations (CATKE closure, WENO advection, SplitRungeKutta3)
#####

tracers = (:T, :S, :dye)
vi = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)
tracer_advection = WENO(order=7, time_discretization=vi, minimum_buffer_upwind_order=1)

shelf_ocean = ocean_simulation(shelf_grid; tracers, tracer_advection,
    boundary_conditions=shelf_boundary_conditions, river_routing=shelf_river_routing, warn=false)

glorys_init_kw = (; dataset=glorys, date=start_date, region=glorys_region)
isnothing(glorys_dir) || (glorys_init_kw = merge(glorys_init_kw, (; dir=glorys_dir)))
set!(shelf_ocean.model, T = Metadatum(:temperature; glorys_init_kw...),
                        S = Metadatum(:salinity;    glorys_init_kw...))

#####
##### Passive dye tracer released from Hinkley Point
#####

@inline hinkley_dye_source(λ, φ, z, t, p) =
    exp(-((λ - p.λ₀)^2 + (φ - p.φ₀)^2) / (2 * p.σ²)) * exp(z / p.depth_scale)

hinkley_dye = Forcing(hinkley_dye_source,
    parameters = (λ₀ = -3.13, φ₀ = 51.21, σ² = 0.005^2, depth_scale = 5.0))

#####
##### EarthSystemModel per level
#####

shelf_model = OceanOnlyModel(shelf_ocean;
    atmosphere=shelf_atmosphere, radiation=shelf_radiation, land=shelf_land)

# OceanOnlyModel calls reconcile_state! → update_state! which refills halos; clamp after
parent(shelf_ocean.model.tracers.S) .= max.(parent(shelf_ocean.model.tracers.S), 0)

bristol_exchanger = StateExchanger(shelf_model)
bristol_sponge    = nested_forcings(bristol_exchanger, bristol_grid;
    rate=1/1hours, width=refinement*width)
bristol_bcs       = nested_boundary_conditions(bristol_exchanger, bristol_grid)

bristol_ocean = ocean_simulation(bristol_grid; tracers, tracer_advection,
    boundary_conditions=bristol_bcs, forcing=bristol_sponge, river_routing=bristol_river_routing, warn=false)
set!(bristol_ocean.model, T = Metadatum(:temperature; glorys_init_kw...),
                          S = Metadatum(:salinity;    glorys_init_kw...))

bristol_model = OceanOnlyModel(bristol_ocean;
    atmosphere=bristol_atmosphere, radiation=bristol_radiation, land=bristol_land)

parent(bristol_ocean.model.tracers.S) .= max.(parent(bristol_ocean.model.tracers.S), 0)

severn_exchanger = StateExchanger(bristol_model)
severn_sponge    = nested_forcings(severn_exchanger, severn_grid;
    rate=1/1hours, width=refinement*width)
severn_bcs       = nested_boundary_conditions(severn_exchanger, severn_grid)

severn_ocean = ocean_simulation(severn_grid; tracers, tracer_advection,
    boundary_conditions=severn_bcs,
    forcing=merge(severn_sponge, (; dye=hinkley_dye)),
    river_routing=severn_river_routing, warn=false)
set!(severn_ocean.model, T = Metadatum(:temperature; glorys_init_kw...),
                          S = Metadatum(:salinity;    glorys_init_kw...))

severn_model = OceanOnlyModel(severn_ocean;
    atmosphere=severn_atmosphere, radiation=severn_radiation,
    land=severn_land)

parent(severn_ocean.model.tracers.S) .= max.(parent(severn_ocean.model.tracers.S), 0)

#####
##### Assemble the three-level nested model
#####

inner_nest = TwoWayNestedModel(bristol_model, severn_model, severn_exchanger; width)
nest       = TwoWayNestedModel(shelf_model, inner_nest, bristol_exchanger; width)

#####
##### Simulation
#####

stop_time = days_to_run * days
simulation = Simulation(nest; Δt, stop_time)

#####
##### Output writers: surface fields every hour, 3D fields every day
#####

for (name, ocean) in pairs((; shelf=shelf_ocean, bristol=bristol_ocean, severn=severn_ocean))
    model = ocean.model
    Nz_model = size(model.grid, 3)

    surface_fields = (; T = model.tracers.T,
                        S = model.tracers.S,
                        u = model.velocities.u,
                        v = model.velocities.v,
                        dye = model.tracers.dye)

    simulation.output_writers[Symbol(name, :_surface)] = JLD2Writer(model, surface_fields;
        filename = joinpath(output_dir, "$(name)_surface"),
        schedule = TimeInterval(1hours),
        indices = (:, :, Nz_model),
        including = [],
        with_halos = false,
        overwrite_files = true)

    simulation.output_writers[Symbol(name, :_eta)] = JLD2Writer(model, (; η = model.free_surface.displacement);
        filename = joinpath(output_dir, "$(name)_eta"),
        schedule = TimeInterval(1hours),
        including = [],
        with_halos = false,
        overwrite_files = true)
#=
    simulation.output_writers[Symbol(name, :_fields)] = JLD2Writer(model, outputs;
        filename = joinpath(output_dir, "$(name)_fields"),
        schedule = TimeInterval(1days),
        with_halos = false,
        overwrite_files = true)=#

end

simulation.output_writers[:checkpoint] = Checkpointer(shelf_ocean.model;
    schedule = TimeInterval(1days),
    dir = output_dir,
    prefix = "checkpoint",
    cleanup = true)

#####
##### Progress
#####

wall_clock = Ref(time_ns())

function progress(sim)
    elapsed = (time_ns() - wall_clock[]) / 1e9
    wall_clock[] = time_ns()
    u = severn_ocean.model.velocities.u
    @printf("  day %7.2f / %.0f  wall: %.1f s  max|u| (severn): %.2f m/s\n",
            sim.model.clock.time / days, days_to_run, elapsed, maximum(abs, u))
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

#####
##### Run
#####

grids = (; shelf=shelf_grid, bristol=bristol_grid, severn=severn_grid)
steps = ceil(Int, stop_time / Δt)
@printf("Severn Estuary: %.0f days, Δt = %.1f s (%d shelf steps)\n", days_to_run, Δt, steps)
for (name, grid) in pairs(grids)
    @printf("  %-12s %3d × %3d × %d, %2.0f%% sea\n",
            name, size(grid)..., 100 * count(sea[name]) / length(sea[name]))
end

run!(simulation)

println("Done. Output in $(output_dir)")
