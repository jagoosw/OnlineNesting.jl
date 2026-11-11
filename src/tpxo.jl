"""
The TPXO10-atlas-v2 barotropic tidal solution on a 1/30° grid.

TPXO is licensed; request the files at https://www.tpxo.net, then unpack the grid and constituent
files into a single directory and pass it as `dir`.
"""
struct TPXO10Atlas
    dir :: String
end

TPXO10Atlas(; dir = joinpath(homedir(), "Downloads", "TPXO10_atlas_v2_nc")) = TPXO10Atlas(dir)

constituent_path(dataset::TPXO10Atlas, variable, constituent) =
    joinpath(dataset.dir, string(variable, "_", lowercase(string(constituent)), "_tpxo10_atlas_30_v2.nc"))

function atlas_window(dataset::TPXO10Atlas, region)
    path = joinpath(dataset.dir, "grid_tpxo10atlas_v2.nc")
    isfile(path) ||
        error("$path not found; TPXO10 must be downloaded by hand from https://www.tpxo.net")

    λᶠ, φᶠ = NCDataset(path) do atlas
        (atlas["lon_u"][:], atlas["lat_v"][:])
    end

    west, east = mod.(region.longitude, 360)
    south, north = region.latitude
    rows = searchsortedlast(φᶠ, south):(searchsortedfirst(φᶠ, north) - 1)

    if west ≤ east
        columns = searchsortedlast(λᶠ, west):(searchsortedfirst(λᶠ, east) - 1)
        longitude = (λᶠ[first(columns)], λᶠ[last(columns) + 1])
    else
        Nλ = length(λᶠ)
        west_columns = searchsortedlast(λᶠ, west):Nλ
        east_columns = 1:(searchsortedfirst(λᶠ, east) - 1)
        columns = (west_columns, east_columns)
        longitude = (λᶠ[first(west_columns)] - 360, λᶠ[last(east_columns) + 1])
    end

    Nc = columns isa Tuple ? length(columns[1]) + length(columns[2]) : length(columns)
    grid = LatitudeLongitudeGrid(size = (Nc, length(rows)),
                                 longitude = longitude,
                                 latitude = (φᶠ[first(rows)], φᶠ[last(rows) + 1]),
                                 topology = (Bounded, Bounded, Flat))

    return (; grid, rows, columns)
end

read_atlas(variable, rows, columns::AbstractRange, face) =
    permutedims(variable[rows, face ? (first(columns):(last(columns) + 1)) : columns])

function read_atlas(variable, rows, columns::Tuple, face)
    west, east = columns
    east_cols = face ? (first(east):(last(east) + 1)) : east
    return permutedims(hcat(variable[rows, west], variable[rows, east_cols]))
end

function atlas_field(window, path, names, (LX, LY), conversion)
    rows = LY === Face ? (first(window.rows):(last(window.rows) + 1)) : window.rows

    real_data, imaginary_data = NCDataset(path) do atlas
        map(name -> read_atlas(atlas[name], rows, window.columns, LX === Face), names)
    end

    mask = Field{LX, LY, Nothing}(window.grid, Bool)
    interior(mask, :, :, 1) .= @. (real_data == 0) & (imaginary_data == 0)

    real_part, imaginary_part = map((real_data, imaginary_data)) do data
        field = Field{LX, LY, Nothing}(window.grid)
        interior(field, :, :, 1) .= conversion .* data
        inpaint_mask!(field, mask)
        field
    end

    constants = Field{LX, LY, Nothing}(window.grid, ComplexF64)
    interior(constants, :, :, 1) .= complex.(interior(real_part, :, :, 1), interior(imaginary_part, :, :, 1))

    return constants
end

"""
$(TYPEDSIGNATURES)

Lowest astronomical tide height [m] relative to mean sea level, computed from the TPXO10 atlas as
the negative sum of tidal elevation amplitudes: LAT ≈ -Σ|Aᵢ|. The result is a `Field` on the
atlas window grid; pass it to [`refer_to_mean_sea_level!`](@ref).
"""
function lowest_astronomical_tide(dataset::TPXO10Atlas, grid;
                                  constituents = keys(earth_tidal_constituents))
    region = (longitude = x_domain(grid) .+ (-1, 1),
              latitude  = y_domain(grid) .+ (-1, 1))
    window = atlas_window(dataset, region)

    lat_height = Field{Center, Center, Nothing}(window.grid)

    for name in constituents
        path = constituent_path(dataset, "h", name)
        η = atlas_field(window, path, ("hRe", "hIm"), (Center, Center), 1e-3)
        interior(lat_height, :, :, 1) .-= abs.(interior(η, :, :, 1))
    end

    return lat_height
end

"""
$(TYPEDSIGNATURES)

Read the harmonic constants of `constituent` at each `(longitude, latitude)` point in `nodes` from
the TPXO10 atlas. Elevations are in meters, eastward and northward transports in meters squared
per second.
"""
function Oceananigans.BoundaryConditions.tidal_atlas_constants(dataset::TPXO10Atlas, nodes, constituent)
    λs = first.(nodes)
    φs = last.(nodes)
    region = (longitude = (minimum(λs) - 1, maximum(λs) + 1), latitude = (minimum(φs) - 1, maximum(φs) + 1))

    window = atlas_window(dataset, region)
    elevation_path = constituent_path(dataset, "h", constituent)
    transport_path = constituent_path(dataset, "u", constituent)

    η_field = atlas_field(window, elevation_path, ("hRe", "hIm"), (Center, Center), 1e-3)
    u_field = atlas_field(window, transport_path, ("uRe", "uIm"), (Face, Center), 1e-4)
    v_field = atlas_field(window, transport_path, ("vRe", "vIm"), (Center, Face), 1e-4)

    λ_west = x_domain(window.grid)[1]
    at_nodes(field) = [interpolate((mod(λ - λ_west, 360) + λ_west, φ), field) for (λ, φ) in nodes]

    return (; sea_surface_height = at_nodes(η_field),
              eastward_transport = at_nodes(u_field),
              northward_transport = at_nodes(v_field))
end
