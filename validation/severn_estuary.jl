# Two-level nested simulation: Bristol Channel shelf → Severn Estuary.
#
# The estuary domain extends east of Gloucester so that all rivers (Severn, Wye,
# Usk, Avon) enter through the interior, not through a lateral boundary.
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
    size = (140, 90, Nz),
    longitude = (-7.0, -2.0),
    latitude = (50.5, 52.5),
    halo = (7, 7, 7),
    topology = (Open, Open, Bounded),
    z)

# The nest extends east of Gloucester (-2.25°W) so the Severn's tidal limit is interior
estuary_grid = nested_grid(shelf_grid;
    longitude = (-5.0, -2.15),
    latitude = (51.1, 52.0),
    depth = 100,
    refinement)

#####
##### Bathymetry (EMODnet referenced to MSL via TPXO lowest astronomical tide)
#####

tpxo = TPXO10Atlas(; dir=tpxo_dir)
lat_field = lowest_astronomical_tide(tpxo, shelf_grid)
lowest_tide(λ, φ) = interpolate((mod(λ, 360), φ), lat_field)

function emodnet_bottom_height(grid; interpolation_passes, major_basins=Inf)
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(grid))
    bottom_height = regrid_bathymetry(grid, metadatum; interpolation_passes, major_basins)
    return refer_to_mean_sea_level!(bottom_height, lowest_tide; minimum_depth)
end

shelf_bottom   = emodnet_bottom_height(shelf_grid;   interpolation_passes=4, major_basins=1)
estuary_bottom = emodnet_bottom_height(estuary_grid; interpolation_passes=2)

blend_parent_bottom!(estuary_bottom, shelf_bottom; width, minimum_depth)
restrict_bottom!(shelf_bottom, estuary_bottom; width, minimum_depth)

sea = map(b -> Array(interior(b, :, :, 1)) .< 0,
          (shelf = shelf_bottom, estuary = estuary_bottom))

immersed(grid, bottom) = ImmersedBoundaryGrid(grid, PartialCellBottom(Array(interior(bottom, :, :, 1))))
shelf_grid, estuary_grid =
    immersed.((shelf_grid, estuary_grid), (shelf_bottom, estuary_bottom))

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

tracers = (:T, :S)

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

estuary_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date,
    region=BoundingBox(estuary_grid), era5_kw...)
estuary_radiation = ERA5PrescribedRadiation(arch; start_date, end_date,
    region=BoundingBox(estuary_grid), era5_kw...)

#####
##### GloFAS river routing on the estuary (optional; set GLOFAS_DIR to enable)
#####

estuary_land = nothing
estuary_river_routing = nothing

if !isnothing(glofas_dir)
    println("Setting up GloFAS river routing on the estuary grid...")
    estuary_land = GloFASPrescribedLand(estuary_grid;
        start_date, end_date, region=BoundingBox(estuary_grid), dir=glofas_dir)
    estuary_river_routing = estuary_land.river_routing
else
    println("Skipping river routing (set GLOFAS_DIR to enable)")
end

#####
##### EarthSystemModel per level (decision D3)
#####

shelf_model = OceanOnlyModel(shelf_ocean;
    atmosphere=shelf_atmosphere, radiation=shelf_radiation)

estuary_exchanger = StateExchanger(shelf_model)
estuary_sponge    = nested_forcings(estuary_exchanger, estuary_grid;
    rate=1/1hour, width=refinement*width)
estuary_bcs       = nested_boundary_conditions(estuary_exchanger, estuary_grid)

estuary_ocean = ocean_simulation(estuary_grid; tracers,
    boundary_conditions=estuary_bcs, forcing=estuary_sponge,
    river_routing=estuary_river_routing, warn=false)
set!(estuary_ocean.model, T = (λ, φ, z) -> initial_temperature(z),
                           S = (λ, φ, z) -> initial_salinity(z))

estuary_model = OceanOnlyModel(estuary_ocean;
    atmosphere=estuary_atmosphere, radiation=estuary_radiation,
    land=estuary_land)

#####
##### Assemble the two-level nested model
#####

nest = TwoWayNestedModel(shelf_model, estuary_model, estuary_exchanger; width)

#####
##### Diagnostics: tide gauge locations
#####

avonmouth     = (-2.71, 51.51)
newport       = (-3.00, 51.55)
mumbles       = (-3.98, 51.57)
ilfracombe    = (-4.11, 51.21)
places = (; avonmouth, newport, mumbles, ilfracombe)

function find_column(grid, sea_mask, (λ, φ))
    λᶠ, φᶠ = Array(λnodes(grid, Face())), Array(φnodes(grid, Face()))
    (λᶠ[1] ≤ λ < λᶠ[end] && φᶠ[1] ≤ φ < φᶠ[end]) || return nothing
    i, j = searchsortedlast(λᶠ, λ), searchsortedlast(φᶠ, φ)
    return sea_mask[i, j] ? (i, j) : nothing
end

grids  = (shelf=shelf_grid, estuary=estuary_grid)
models = (shelf=shelf_ocean.model, estuary=estuary_ocean.model)

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

Δt = shelf_ocean.Δt
stop_time = days_to_run * days
steps = ceil(Int, stop_time / Δt)
record_interval = max(1, steps ÷ 1000)

println("Running the Severn Estuary case for $(days_to_run) days, Δt = $(Δt) s ($(steps) shelf steps)")
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
    t = series.time
    for place in keys(places)
        println("  $(place):")
        η_shelf   = series.η.shelf[place]
        η_estuary = series.η.estuary[place]
        if !any(isnan, η_estuary) && !any(isnan, η_shelf)
            rms = sqrt(sum((η_shelf .- η_estuary).^2) / length(η_shelf))
            @printf("    η shelf vs estuary: range %.2f vs %.2f m, rms difference %.4f m\n",
                    maximum(η_shelf) - minimum(η_shelf),
                    maximum(η_estuary) - minimum(η_estuary), rms)
        end
        s_shelf   = series.speed.shelf[place]
        s_estuary = series.speed.estuary[place]
        if !any(isnan, s_estuary) && !any(isnan, s_shelf)
            @printf("    speed shelf vs estuary: peak %.2f vs %.2f m/s\n",
                    maximum(s_shelf), maximum(s_estuary))
        end
    end
end

println("\nAgreement between levels:")
report_agreement(series, places)

println("\nDone.")
