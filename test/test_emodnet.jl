using Test
using Oceananigans
using Oceananigans.Grids: λnodes, φnodes
using NumericalEarth
using NumericalEarth.DataWrangling: BoundingBox, Metadatum, metadata_path, native_grid
using NCDatasets
using OnlineNesting

# Skomer Island and the sea around it
longitude, latitude = (-5.33, -5.26), (51.72, 51.75)

@testset "EMODnet2024" begin
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), dir=mktempdir(),
                          region=BoundingBox(; longitude, latitude))
    grid = LatitudeLongitudeGrid(size=(28, 12, 1); longitude, latitude, z=(-100, 0))
    bottom_height = Array(interior(regrid_bathymetry(grid, metadatum; cache=false, major_basins=Inf), :, :, 1))

    NCDataset(metadata_path(metadatum)) do dataset
        source_grid = native_grid(metadatum)
        @test maximum(abs, dataset["longitude"][:] .- λnodes(source_grid, Center())) < 1 / 9600
        @test maximum(abs, dataset["latitude"][:] .- φnodes(source_grid, Center())) < 1 / 9600
        @test size(dataset["elevation"]) == size(source_grid)[1:2]
    end

    @test -60 < minimum(bottom_height) < -30
    @test 0.1 < count(>(0), bottom_height) / length(bottom_height) < 0.5
end
