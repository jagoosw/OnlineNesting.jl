"""
    EMODnet2024

The EMODnet 2024 digital terrain model of the European seas: elevation [m] relative to mean sea
level, 1/16 arc-minute (≈ 115 m) apart, from 36°W to 43°E and 15°N to 90°N. Land without data
is given a height of 10 m.

The data comes in tiles of about 100 to 200 MB, which are downloaded whole and kept. A
`Metadatum` of this dataset takes a `BoundingBox` as its `region` and reads that window of them.

Data source: https://emodnet.ec.europa.eu/en/bathymetry
"""
struct EMODnet2024 <: AbstractStaticBathymetry end

const EMODnetMetadatum = Metadatum{<:EMODnet2024}

const EMODnet_url = "https://downloads.emodnet-bathymetry.eu/v12"

DataWrangling.default_download_directory(::EMODnet2024) = download_cache("EMODnet")
DataWrangling.longitude_interfaces(::EMODnet2024) = (-36, 43)
DataWrangling.latitude_interfaces(::EMODnet2024) = (15, 90)
Base.size(::EMODnet2024) = (75840, 72000, 1)

DataWrangling.dataset_variable_name(::EMODnetMetadatum) = "elevation"
DataWrangling.metadata_filename(::EMODnet2024, name, date, region) = "EMODnet2024_mean_sea_level_" * bounding_box_suffix(region) * ".nc"

# Copy from an ESRI ASCII tile the part within the window that starts `west` cells east of 36°W
# and `north` cells south of 90°N
function read_tile!(elevation, tile_path, west, north)
    Nx, Ny = size(elevation)

    open(tile_path) do tile
        columns, rows, west_edge, south_edge, _, no_data = [parse(Float64, split(readline(tile))[2]) for _ in 1:6]
        tile_west = round(Int, (west_edge + 36) * 960)
        tile_north = round(Int, (90 - south_edge) * 960) - Int(rows)
        i = max(1, tile_west - west + 1):min(Nx, tile_west - west + Int(columns))

        for (row, line) in enumerate(eachline(tile))
            j = Ny - (tile_north + row - 1 - north)
            j > Ny && continue
            j < 1 && break
            heights = parse.(Float32, view(split(line), i .+ (west - tile_west)))
            elevation[i, j] .= ifelse.(heights .== no_data, 10, heights)
        end
    end

    return nothing
end

function Downloads.download(metadatum::EMODnetMetadatum)
    path = metadata_path(metadatum)
    isfile(path) && return path

    grid = native_grid(metadatum)
    Nx, Ny, _ = size(grid)
    west = round(Int, mod(x_domain(grid)[1] + 36, 360) * 960)
    north = round(Int, (90 - y_domain(grid)[2]) * 960)
    elevation = zeros(Float32, Nx, Ny)

    # Tiles are 9480 × 9000 cells, lettered southward from 90°N and numbered eastward from 36°W
    for row in north ÷ 9000:(north + Ny - 1) ÷ 9000, column in west ÷ 9480:(west + Nx - 1) ÷ 9480
        tile = string('A' + row, column + 1)
        archive = joinpath(dirname(path), "$(tile)_2024.msl.zip")
        isfile(archive) || download_with_retries("$EMODnet_url/$(tile)_2024.msl.zip", archive)

        directory = unzip(archive, mktempdir())
        read_tile!(elevation, joinpath(directory, "$(tile)_MSL_2024.asc"), west, north)
        rm(directory, recursive=true)
    end

    staging_path = tempname() * ".nc"

    NCDataset(staging_path, "c") do dataset
        defVar(dataset, "longitude", collect(λnodes(grid, Center())), ("longitude",))
        defVar(dataset, "latitude", collect(φnodes(grid, Center())), ("latitude",))
        defVar(dataset, "elevation", elevation, ("longitude", "latitude"))
    end

    return mv(staging_path, path)
end
