"""
    EMODnet2024

The EMODnet 2024 digital terrain model of the European seas: elevation [m] relative to the lowest
astronomical tide, 1/16 arc-minute (≈ 115 m) apart, from 36°W to 43°E and 15°N to 90°N. Land
without data is given a height of 10 m.

The data is read in windows from the EMODnet ERDDAP server, so a `Metadatum` of this dataset
takes a `BoundingBox` as its `region`.

Data source: https://emodnet.ec.europa.eu/en/bathymetry
"""
struct EMODnet2024 <: AbstractStaticBathymetry end

const EMODnetMetadatum = Metadatum{<:EMODnet2024}

const EMODnet_url = "https://erddap.emodnet.eu/erddap/griddap/bathymetry_dtm_2024.nc"

DataWrangling.default_download_directory(::EMODnet2024) = download_cache("EMODnet")
DataWrangling.longitude_interfaces(::EMODnet2024) = (-36, 43)
DataWrangling.latitude_interfaces(::EMODnet2024) = (15, 90)
Base.size(::EMODnet2024) = (75840, 72000, 1)

DataWrangling.dataset_variable_name(::EMODnetMetadatum) = "elevation"
DataWrangling.metadata_filename(::EMODnet2024, name, date, region) = "EMODnet2024_" * bounding_box_suffix(region) * ".nc"

function Downloads.download(metadatum::EMODnetMetadatum)
    path = metadata_path(metadatum)
    isfile(path) && return path

    grid = native_grid(metadatum)
    Nx, Ny, _ = size(grid)
    i = round(Int, mod(x_domain(grid)[1] + 36, 360) * 960)
    j = round(Int, (y_domain(grid)[1] - 15) * 960)

    # ERDDAP indices start at zero, with latitude first
    staging_path = download_with_retries("$EMODnet_url?elevation[$j:1:$(j + Ny - 1)][$i:1:$(i + Nx - 1)]", tempname() * ".nc")

    NCDataset(staging_path, "a") do dataset
        elevation = dataset["elevation"]
        elevation[:, :] = coalesce.(elevation[:, :], 10)
    end

    return mv(staging_path, path)
end
