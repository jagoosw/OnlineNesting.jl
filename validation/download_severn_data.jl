# Download GLORYS and GloFAS data for the Severn Estuary three-level nesting case (CPU only).
#
#     julia --project=test validation/download_severn_data.jl
#
# Already-downloaded files are skipped, so it can be rerun after an interruption.
# ERA5 and EMODnet bathymetry are downloaded automatically when the simulation runs.
# TPXO10 is licensed and must be obtained separately (https://www.tpxo.net).
#
# Copernicus Marine credentials (GLORYS): run `copernicusmarine login` once, or set
#   COPERNICUSMARINE_SERVICE_USERNAME and COPERNICUSMARINE_SERVICE_PASSWORD.
#
# CDS credentials (GloFAS): create ~/.cdsapirc with your CDS API key from
#   https://cds.climate.copernicus.eu, and accept the GloFAS licence there.

data_directory = get(ENV, "NUMERICALEARTH_DATA_DIRECTORY", "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data")
ENV["NUMERICALEARTH_DATA_DIRECTORY"] = data_directory

using NumericalEarth
using CopernicusMarine
using CDSAPI
using Downloads
using Dates: DateTime, Day

#####
##### Configuration — must match severn_estuary.jl
#####

start_date = DateTime(2018, 1, 1)
end_date   = DateTime(2021, 1, 1)

# Shelf domain with margin for GLORYS interpolation
shelf_longitude = (-12.0, 5.0)
shelf_latitude  = (46.0, 62.0)
glorys_region = BoundingBox(longitude = (shelf_longitude[1] - 0.5, shelf_longitude[2] + 0.5),
                            latitude  = (shelf_latitude[1]  - 0.5, shelf_latitude[2]  + 0.5))

# Severn domain for GloFAS river routing (covers Severn, Wye, Usk, Avon mouths)
severn_longitude = (-3.5, -2.5)
severn_latitude  = (51.1, 51.75)
glofas_region = BoundingBox(; longitude=severn_longitude, latitude=severn_latitude)

#####
##### GLORYS daily reanalysis (T, S, u, v, η)
#####

glorys = GLORYSDaily()
glorys_dates = start_date:Day(1):(end_date - Day(1))

@sync begin
    for variable in (:temperature, :salinity, :u_velocity, :v_velocity, :free_surface)
        @async begin
            @info "GLORYS $variable"
            Downloads.download(Metadata(variable; dataset=glorys, dates=glorys_dates, region=glorys_region))
            @info "GLORYS $variable done"
        end
    end

    @async begin
        @info "GloFAS river discharge"
        Downloads.download(Metadata(:river_discharge; dataset=GloFASReanalysis(),
                                    start_date, end_date, region=glofas_region))
        @info "GloFAS done"
    end
end

@info "Done"
