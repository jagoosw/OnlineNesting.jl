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
using Downloads: Downloads
using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: TracerReservoir, tidal_boundary_conditions
using Oceananigans.Fields: interpolate
using Oceananigans.Grids: Open, MutableVerticalDiscretization, λnodes, φnodes
using NumericalEarth
using NumericalEarth.DataWrangling: BoundingBox, GLORYSDaily, Metadata, MetadataSet
using NumericalEarth.NestedModels: Interpolated
using OnlineNesting
using Printf

#####
##### Configuration
#####

arch = CPU()
days_to_run = isempty(ARGS) ? 30 : parse(Float64, ARGS[1])
start_date = DateTime(2020, 1, 1)
end_date   = DateTime(2022, 1, 1)
dates = (start_date, end_date)

tpxo_dir   = get(ENV, "TPXO_DIR",   joinpath(homedir(), "Downloads", "TPXO10_atlas_v2_nc"))
glorys_dir = get(ENV, "GLORYS_DIR", nothing)
era5_dir   = get(ENV, "ERA5_DIR",   nothing)
glofas_dir = get(ENV, "GLOFAS_DIR", nothing)

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
Δt = result.timesteps[1]

#####
##### Bathymetry (EMODnet referenced to MSL via TPXO lowest astronomical tide)
#####

tpxo = TPXO10Atlas(; dir=tpxo_dir)
# LAT computed on bristol_grid: the shelf grid spans 0° longitude,
# triggering the atlas_window wrapping bug in TPXO interpolation
lat_field = lowest_astronomical_tide(tpxo, bristol_grid)
lowest_tide(λ, φ) = interpolate((mod(λ, 360), φ), lat_field)

function emodnet_bottom_height(grid; interpolation_passes, major_basins=Inf, lat_correction=lowest_tide)
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(grid))
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
shelf_region = BoundingBox(shelf_grid)
glorys_kw = (; dataset=glorys, start_date, end_date, region=shelf_region)
isnothing(glorys_dir) || (glorys_kw = merge(glorys_kw, (; dir=glorys_dir)))
fts_kw = (; time_indices_in_memory=30)

mset = MetadataSet(:temperature, :salinity, :u_velocity, :v_velocity; glorys_kw...)
Downloads.download(mset)

glorys_T = FieldTimeSeries(Metadata(:temperature; glorys_kw...), arch; fts_kw...)
glorys_S = FieldTimeSeries(Metadata(:salinity;    glorys_kw...), arch; fts_kw...)
glorys_u = FieldTimeSeries(Metadata(:u_velocity;  glorys_kw...), arch; fts_kw...)
glorys_v = FieldTimeSeries(Metadata(:v_velocity;  glorys_kw...), arch; fts_kw...)

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
##### Ocean simulations (CATKE closure, WENO advection, SplitRungeKutta3)
#####

tracers = (:T, :S, :dye)

shelf_ocean = ocean_simulation(shelf_grid; tracers,
    boundary_conditions=shelf_boundary_conditions, warn=false)

initial_temperature(z) = 10 + 8 * exp(z / 50)
initial_salinity(z) = 35.0 + 0.1 * exp(z / 100)
set!(shelf_ocean.model, T = (λ, φ, z) -> initial_temperature(z),
                         S = (λ, φ, z) -> initial_salinity(z))

#####
##### ERA5 atmosphere and radiation per level
#####

era5_kw = isnothing(era5_dir) ? (;) : (; dir=era5_dir)

shelf_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date,
    region=BoundingBox(shelf_grid), era5_kw...)
shelf_radiation = ERA5PrescribedRadiation(arch; start_date, end_date,
    region=BoundingBox(shelf_grid), era5_kw...)

bristol_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date,
    region=BoundingBox(bristol_grid), era5_kw...)
bristol_radiation = ERA5PrescribedRadiation(arch; start_date, end_date,
    region=BoundingBox(bristol_grid), era5_kw...)

severn_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date,
    region=BoundingBox(severn_grid), era5_kw...)
severn_radiation = ERA5PrescribedRadiation(arch; start_date, end_date,
    region=BoundingBox(severn_grid), era5_kw...)

#####
##### GloFAS river routing on the Severn (optional; set GLOFAS_DIR to enable)
#####

severn_land = nothing
severn_river_routing = nothing

if !isnothing(glofas_dir)
    println("Setting up GloFAS river routing on the Severn grid...")
    severn_land = GloFASPrescribedLand(severn_grid;
        start_date, end_date, region=BoundingBox(severn_grid), dir=glofas_dir)
    severn_river_routing = severn_land.river_routing
else
    println("Skipping river routing (set GLOFAS_DIR to enable)")
end

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
    atmosphere=shelf_atmosphere, radiation=shelf_radiation)

bristol_exchanger = StateExchanger(shelf_model)
bristol_sponge    = nested_forcings(bristol_exchanger, bristol_grid;
    rate=1/1hour, width=refinement*width)
bristol_bcs       = nested_boundary_conditions(bristol_exchanger, bristol_grid)

bristol_ocean = ocean_simulation(bristol_grid; tracers,
    boundary_conditions=bristol_bcs, forcing=bristol_sponge, warn=false)
set!(bristol_ocean.model, T = (λ, φ, z) -> initial_temperature(z),
                           S = (λ, φ, z) -> initial_salinity(z))

bristol_model = OceanOnlyModel(bristol_ocean;
    atmosphere=bristol_atmosphere, radiation=bristol_radiation)

severn_exchanger = StateExchanger(bristol_model)
severn_sponge    = nested_forcings(severn_exchanger, severn_grid;
    rate=1/1hour, width=refinement*width)
severn_bcs       = nested_boundary_conditions(severn_exchanger, severn_grid)

severn_ocean = ocean_simulation(severn_grid; tracers,
    boundary_conditions=severn_bcs,
    forcing=merge(severn_sponge, (; dye=hinkley_dye)),
    river_routing=severn_river_routing, warn=false)
set!(severn_ocean.model, T = (λ, φ, z) -> initial_temperature(z),
                          S = (λ, φ, z) -> initial_salinity(z))

severn_model = OceanOnlyModel(severn_ocean;
    atmosphere=severn_atmosphere, radiation=severn_radiation,
    land=severn_land)

#####
##### Assemble the three-level nested model
#####

inner_nest = TwoWayNestedModel(bristol_model, severn_model, severn_exchanger; width)
nest       = TwoWayNestedModel(shelf_model, inner_nest, bristol_exchanger; width)

#####
##### Diagnostics: tide gauge locations
#####

avonmouth     = (-2.71, 51.51)
newport       = (-3.00, 51.55)
mumbles       = (-3.98, 51.57)
ilfracombe    = (-4.11, 51.21)
hinkley_point = (-3.13, 51.21)
places = (; avonmouth, newport, mumbles, ilfracombe, hinkley_point)

function find_column(grid, sea_mask, (λ, φ))
    λᶠ, φᶠ = Array(λnodes(grid, Face())), Array(φnodes(grid, Face()))
    (λᶠ[1] ≤ λ < λᶠ[end] && φᶠ[1] ≤ φ < φᶠ[end]) || return nothing
    i, j = searchsortedlast(λᶠ, λ), searchsortedlast(φᶠ, φ)
    return sea_mask[i, j] ? (i, j) : nothing
end

grids  = (shelf=shelf_grid, bristol=bristol_grid, severn=severn_grid)
models = (shelf=shelf_ocean.model, bristol=bristol_ocean.model, severn=severn_ocean.model)

columns = map((g, s) -> map(p -> find_column(g, s, p), places), values(grids), values(sea))
columns = NamedTuple{keys(grids)}(columns)

series = (time  = Float64[],
          η     = map(_ -> map(_ -> Float64[], places), grids),
          speed = map(_ -> map(_ -> Float64[], places), grids))

function record!(series, models, columns, time)
    push!(series.time, time)
    for name in keys(models)
        model = models[name]
        η_data = Array(interior(model.free_surface.displacement, :, :, 1))
        Nz_model = size(model.grid, 3)
        u_data = Array(interior(model.velocities.u, :, :, Nz_model))
        v_data = Array(interior(model.velocities.v, :, :, Nz_model))
        speed_data = sqrt.(((u_data[1:end-1, :] .+ u_data[2:end, :]) ./ 2) .^ 2 .+
                           ((v_data[:, 1:end-1] .+ v_data[:, 2:end]) ./ 2) .^ 2)
        for place in keys(places)
            col = columns[name][place]
            push!(series.η[name][place],     isnothing(col) ? NaN : η_data[col...])
            push!(series.speed[name][place], isnothing(col) ? NaN : speed_data[col...])
        end
    end
end

#####
##### Run
#####

stop_time = days_to_run * days
steps = ceil(Int, stop_time / Δt)
record_interval = max(1, steps ÷ 1000)

println("Running the Severn Estuary case for $(days_to_run) days, Δt = $(round(Δt, digits=1)) s ($(steps) shelf steps)")
for name in keys(grids)
    g = grids[name]
    idx = findfirst(==(name), keys(grids))
    @printf("  %-12s %3d × %3d × %d, %2.0f%% sea\n",
            name, size(g)..., 100 * count(values(sea)[idx]) / length(values(sea)[idx]))
end

for step in 1:steps
    time_step!(nest, Δt)

    if step % record_interval == 0
        record!(series, models, columns, nest.parent.clock.time)
        t_days = nest.parent.clock.time / days
        @printf("\r  day %7.1f / %7.1f", t_days, days_to_run)
    end
end
println()

#####
##### Report agreement between levels
#####

function report_agreement(series, places)
    for place in keys(places)
        println("  $(place):")
        for (parent_name, child_name) in ((:shelf, :bristol), (:bristol, :severn))
            η_parent = series.η[parent_name][place]
            η_child  = series.η[child_name][place]
            if !any(isnan, η_child) && !any(isnan, η_parent)
                rms = sqrt(sum((η_parent .- η_child).^2) / length(η_parent))
                @printf("    η %s vs %s: range %.2f vs %.2f m, rms difference %.4f m\n",
                        parent_name, child_name,
                        maximum(η_parent) - minimum(η_parent),
                        maximum(η_child) - minimum(η_child), rms)
            end
            s_parent = series.speed[parent_name][place]
            s_child  = series.speed[child_name][place]
            if !any(isnan, s_child) && !any(isnan, s_parent)
                @printf("    speed %s vs %s: peak %.2f vs %.2f m/s\n",
                        parent_name, child_name,
                        maximum(s_parent), maximum(s_child))
            end
        end
    end
end

println("\nAgreement between levels:")
report_agreement(series, places)

println("\nDone.")
