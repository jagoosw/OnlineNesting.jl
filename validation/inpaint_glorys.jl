# Pre-inpaint all GLORYS daily fields so the simulation finds cached .jld2 files
# instead of inpainting on the fly.
#
#     julia --project=test -t auto validation/inpaint_glorys.jl
#
# Expects GLORYS data already downloaded by download_severn_data.jl.

using Dates
using CopernicusMarine
using Oceananigans
using NumericalEarth
using NumericalEarth.DataWrangling: GLORYSDaily, Metadatum, BoundingBox, default_inpainting
using Printf

glorys_dir = get(ENV, "GLORYS_DIR",
    "/nfs/roberts/pi/pi_ey239/js5256/OAEMIP/data/GLORYS/")

arch = CPU()
glorys = GLORYSDaily()
region = BoundingBox(longitude = (-12.5, 5.5), latitude = (45.5, 62.5))

start_date = DateTime(2018, 1, 1)
end_date   = DateTime(2021, 1, 1)
all_dates  = start_date:Day(1):(end_date - Day(1))

variables = (:temperature, :salinity, :u_velocity, :v_velocity)

for variable in variables
    @info "Inpainting $variable ($(length(all_dates)) dates)..."
    t₀ = time_ns()
    for (n, date) in enumerate(all_dates)
        metadatum = Metadatum(variable; dataset=glorys, date, region, dir=glorys_dir)
        Field(metadatum, arch; cache_inpainted_data=true)
        n % 50 == 0 && @printf("  %s  %d / %d\n", variable, n, length(all_dates))
    end
    elapsed = (time_ns() - t₀) / 1e9
    @printf("  %s done in %.0f s\n", variable, elapsed)
end

@info "All inpainting cached."
