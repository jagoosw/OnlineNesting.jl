# Three-level nested ocean simulation of the NW European shelf, Pembrokeshire approaches
# and St Brides Bay, with realistic forcing (ERA5 atmosphere, GLORYS reanalysis open boundary,
# TPXO10 tides) and CATKE turbulence.
#
# Prerequisites:
#   - TPXO10 atlas files in TPXO_DIR (licensed; request at https://www.tpxo.net)
#   - Copernicus Marine credentials: COPERNICUSMARINE_SERVICE_USERNAME, _PASSWORD
#   - ERA5 data auto-downloads; set ERA5_DIR for a shared cache
#
#   julia --project=test -t auto validation/pembrokeshire.jl [days, default 730]
#
# For GPU, load CUDA before this script and set arch = GPU() below.

using Dates
using Downloads: Downloads
using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: TracerReservoir, tidal_boundary_conditions
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
days_to_run = isempty(ARGS) ? 730 : parse(Float64, ARGS[1])
start_date = DateTime(2020, 1, 1)
end_date = DateTime(2022, 1, 1)
dates = (start_date, end_date)

tpxo_dir = get(ENV, "TPXO_DIR", joinpath(homedir(), "Downloads", "TPXO10_atlas_v2_nc"))
glorys_dir = get(ENV, "GLORYS_DIR", nothing)
era5_dir = get(ENV, "ERA5_DIR", nothing)

refinement = 5
width = 3
minimum_depth = 10
Nz = 50

#####
##### Vertical grid: 50 surface-intensified layers to 250 m
#####

z_faces = [-250 * (1 - k / Nz)^1.8 for k in 0:Nz]
z = MutableVerticalDiscretization(z_faces)

#####
##### Horizontal grids
#####

shelf_grid = LatitudeLongitudeGrid(arch;
    size = (220, 245, Nz),
    longitude = (-10.5, -2.5),
    latitude = (49.0, 54.5),
    halo = (7, 7, 7),
    topology = (Open, Open, Bounded),
    z)

approaches_grid = nested_grid(shelf_grid; longitude=(-6.4, -4.6), latitude=(51.3, 52.4), depth=150, refinement)
bay_grid = nested_grid(approaches_grid; longitude=(-5.50, -5.06), latitude=(51.68, 51.94), depth=95, refinement)

#####
##### Bathymetry
#####

function emodnet_bottom_height(grid; interpolation_passes, major_basins=Inf)
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(grid))
    bottom_height = regrid_bathymetry(grid, metadatum; interpolation_passes, major_basins)
    return refer_to_mean_sea_level!(bottom_height, 0; minimum_depth)
end

shelf_bottom = emodnet_bottom_height(shelf_grid; interpolation_passes=4, major_basins=1)
approaches_bottom = emodnet_bottom_height(approaches_grid; interpolation_passes=2)
bay_bottom = emodnet_bottom_height(bay_grid; interpolation_passes=1)

blend_parent_bottom!(approaches_bottom, shelf_bottom; width, minimum_depth)
blend_parent_bottom!(bay_bottom, approaches_bottom; width, minimum_depth)
restrict_bottom!(approaches_bottom, bay_bottom; width, minimum_depth)
restrict_bottom!(shelf_bottom, approaches_bottom; width, minimum_depth)

sea = map(b -> Array(interior(b, :, :, 1)) .< 0,
          (shelf = shelf_bottom, approaches = approaches_bottom, bay = bay_bottom))

immersed(grid, bottom) = ImmersedBoundaryGrid(grid, PartialCellBottom(Array(interior(bottom, :, :, 1))))
shelf_grid, approaches_grid, bay_grid =
    immersed.((shelf_grid, approaches_grid, bay_grid),
              (shelf_bottom, approaches_bottom, bay_bottom))

#####
##### Tidal forcing on the shelf boundary
#####

harmonics = earth_tidal_harmonics(start_date; ramp_time=6hours)
tidal_bcs = tidal_boundary_conditions(shelf_grid, harmonics, TPXO10Atlas(; dir=tpxo_dir))

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
glorys_S = FieldTimeSeries(Metadata(:salinity; glorys_kw...), arch; fts_kw...)
glorys_u = FieldTimeSeries(Metadata(:u_velocity; glorys_kw...), arch; fts_kw...)
glorys_v = FieldTimeSeries(Metadata(:v_velocity; glorys_kw...), arch; fts_kw...)

reservoir(source) = ValueBoundaryCondition(Interpolated(source); scheme=TracerReservoir())
normal(source) = NormalFlowBoundaryCondition(Interpolated(source))
zero_reservoir = ValueBoundaryCondition(0; scheme=TracerReservoir())

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
    P = FieldBoundaryConditions(west=zero_reservoir, east=zero_reservoir,
                                south=zero_reservoir, north=zero_reservoir),
)

#####
##### Ocean simulations (CATKE closure, WENO advection, SplitRungeKutta3)
#####

tracers = (:T, :S, :P)

shelf_ocean = ocean_simulation(shelf_grid; tracers, boundary_conditions=shelf_boundary_conditions,
                               warn=false)

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

approaches_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date,
    region=BoundingBox(approaches_grid), era5_kw...)
approaches_radiation = ERA5PrescribedRadiation(arch; start_date, end_date,
    region=BoundingBox(approaches_grid), era5_kw...)

bay_atmosphere = ERA5PrescribedAtmosphere(arch; start_date, end_date,
    region=BoundingBox(bay_grid), era5_kw...)
bay_radiation = ERA5PrescribedRadiation(arch; start_date, end_date,
    region=BoundingBox(bay_grid), era5_kw...)

#####
##### EarthSystemModel per level (decision D3)
#####

shelf_model = OceanOnlyModel(shelf_ocean; atmosphere=shelf_atmosphere, radiation=shelf_radiation)

# Approaches: nested from shelf
approaches_exchanger = StateExchanger(shelf_model)
approaches_sponge = nested_forcings(approaches_exchanger, approaches_grid; rate=1/1hour, width=refinement*width)
approaches_bcs = nested_boundary_conditions(approaches_exchanger, approaches_grid)

approaches_ocean = ocean_simulation(approaches_grid; tracers,
    boundary_conditions=approaches_bcs, forcing=approaches_sponge, warn=false)
set!(approaches_ocean.model, T = (λ, φ, z) -> initial_temperature(z),
                              S = (λ, φ, z) -> initial_salinity(z))

approaches_model = OceanOnlyModel(approaches_ocean;
    atmosphere=approaches_atmosphere, radiation=approaches_radiation)

# Bay: nested from approaches
bay_exchanger = StateExchanger(approaches_model)
bay_sponge = nested_forcings(bay_exchanger, bay_grid; rate=1/1hour, width=refinement*width)
bay_bcs = nested_boundary_conditions(bay_exchanger, bay_grid)

bay_ocean = ocean_simulation(bay_grid; tracers,
    boundary_conditions=bay_bcs, forcing=bay_sponge, warn=false)
set!(bay_ocean.model, T = (λ, φ, z) -> initial_temperature(z),
                       S = (λ, φ, z) -> initial_salinity(z))

bay_model = OceanOnlyModel(bay_ocean; atmosphere=bay_atmosphere, radiation=bay_radiation)

#####
##### Passive tracer released in the bay
#####

P_bay = bay_ocean.model.tracers.P
λc = Array(λnodes(bay_grid, Center()))
φc = Array(φnodes(bay_grid, Center()))
release_mask = @. (λc > -5.35) & (λc < -5.20) & (φc' > 51.75) & (φc' < 51.85)
interior(P_bay)[:, :, end] .= release_mask

#####
##### Assemble the three-level nested model
#####

inner_nest = TwoWayNestedModel(approaches_model, bay_model, bay_exchanger; width)
nest = TwoWayNestedModel(shelf_model, inner_nest, approaches_exchanger; width)

#####
##### Diagnostics
#####

st_brides_bay = (-5.25, 51.80)
ramsey_sound = (-5.3225, 51.868)
places = (; st_brides_bay, ramsey_sound)

function find_column(grid, sea_mask, (λ, φ))
    λᶠ, φᶠ = Array(λnodes(grid, Face())), Array(φnodes(grid, Face()))
    (λᶠ[1] ≤ λ < λᶠ[end] && φᶠ[1] ≤ φ < φᶠ[end]) || return nothing
    i, j = searchsortedlast(λᶠ, λ), searchsortedlast(φᶠ, φ)
    return sea_mask[i, j] ? (i, j) : nothing
end

grids = (shelf=shelf_grid, approaches=approaches_grid, bay=bay_grid)
models = (shelf=shelf_ocean.model, approaches=approaches_ocean.model, bay=bay_ocean.model)

columns = map((g, s) -> map(p -> find_column(g, s, p), places), values(grids), values(sea))
columns = NamedTuple{keys(grids)}(columns)

series = (time = Float64[],
          η = map(_ -> map(_ -> Float64[], places), grids),
          speed = map(_ -> map(_ -> Float64[], places), grids))

function record!(series, models, columns, time)
    push!(series.time, time)
    for name in keys(models)
        model = models[name]
        η_data = Array(interior(model.free_surface.displacement, :, :, 1))
        Nz = size(model.grid, 3)
        u_data = Array(interior(model.velocities.u, :, :, Nz))
        v_data = Array(interior(model.velocities.v, :, :, Nz))
        speed_data = sqrt.(((u_data[1:end-1, :] .+ u_data[2:end, :]) ./ 2) .^ 2 .+
                           ((v_data[:, 1:end-1] .+ v_data[:, 2:end]) ./ 2) .^ 2)
        for place in keys(places)
            col = columns[name][place]
            push!(series.η[name][place], isnothing(col) ? NaN : η_data[col...])
            push!(series.speed[name][place], isnothing(col) ? NaN : speed_data[col...])
        end
    end
end

#####
##### Run with feedback on
#####

Δt = shelf_ocean.Δt
stop_time = days_to_run * days
steps = ceil(Int, stop_time / Δt)
record_interval = max(1, steps ÷ 1000)

println("Running the Pembrokeshire case for $(days_to_run) days, Δt = $(Δt) s ($(steps) shelf steps)")
for (name, g, s) in zip(keys(grids), values(grids), values(sea))
    @printf("  %-12s %3d × %3d × %d, %2.0f%% sea\n", name, size(g)..., 100 * count(s) / length(s))
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
        for (outer, inner) in [(:shelf, :approaches), (:approaches, :bay)]
            η_out = series.η[outer][place]
            η_in = series.η[inner][place]
            if !any(isnan, η_in) && !any(isnan, η_out)
                rms = sqrt(sum((η_out .- η_in).^2) / length(η_out))
                @printf("    η %s vs %s: range %.2f vs %.2f m, rms difference %.4f m\n",
                        outer, inner, maximum(η_out) - minimum(η_out),
                        maximum(η_in) - minimum(η_in), rms)
            end
            s_out = series.speed[outer][place]
            s_in = series.speed[inner][place]
            if !any(isnan, s_in) && !any(isnan, s_out)
                @printf("    speed %s vs %s: peak %.2f vs %.2f m/s\n",
                        outer, inner, maximum(s_out), maximum(s_in))
            end
        end
    end
end

println("\nAgreement between levels (feedback on):")
report_agreement(series, places)

#####
##### Check: passive tracer appears on parent grids
#####

for (name, model) in pairs(models)
    P = model.tracers.P
    max_P = maximum(P)
    total_P = sum(interior(P))
    @printf("Passive tracer %s: max = %.4e, total = %.4e\n", name, max_P, total_P)
end

# Composite conservation is exercised by the idealized tests (milestone 7) and verified to
# round-off on z-star. For the feedback-off comparison, rerun with `width = 0` in the
# TwoWayNestedModel constructors.

println("\nDone.")
