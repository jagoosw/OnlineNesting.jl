# EMODnet bathymetry on the Bristol Channel shelf and Severn Estuary nest,
# reconciled between levels. Draws each level's depth with the outline of the
# nest, and reports column volume agreement.

using Oceananigans
using Oceananigans.Fields: interpolate
using Oceananigans.Grids: static_column_depthᶜᶜᵃ, x_domain, y_domain, λnodes, φnodes
using Oceananigans.Operators: Azᶜᶜᶜ
using NumericalEarth
using NumericalEarth.DataWrangling: BoundingBox, Metadatum
using OnlineNesting
using CairoMakie
using Printf

refinement    = 5
width         = 3
minimum_depth = 10
Nz            = 50
z_faces = [-250 * (1 - k / Nz)^1.8 for k in 0:Nz]

shelf_grid = LatitudeLongitudeGrid(size=(140, 90, Nz),
    longitude=(-7.0, -2.0), latitude=(50.5, 52.5), halo=(7, 7, 4); z=z_faces)

estuary_grid = nested_grid(shelf_grid;
    longitude=(-5.0, -2.15), latitude=(51.1, 52.0), depth=100, refinement)

data_dir = joinpath(dirname(@__DIR__), "data")
tpxo_dir = get(ENV, "TPXO_DIR", joinpath(homedir(), "Downloads", "TPXO10_atlas_v2_nc"))

tpxo = TPXO10Atlas(; dir=tpxo_dir)
lat_field = lowest_astronomical_tide(tpxo, shelf_grid)
lowest_tide(λ, φ) = interpolate((mod(λ, 360), φ), lat_field)

function emodnet_bottom_height(grid; interpolation_passes, major_basins=Inf)
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(grid), dir=data_dir)
    bottom_height = regrid_bathymetry(grid, metadatum; minimum_depth, interpolation_passes, major_basins)
    return refer_to_mean_sea_level!(bottom_height, lowest_tide; minimum_depth)
end

println("Downloading and regridding bathymetry...")
shelf_bottom   = emodnet_bottom_height(shelf_grid;   interpolation_passes=4, major_basins=1)
estuary_bottom = emodnet_bottom_height(estuary_grid; interpolation_passes=2)

blend_parent_bottom!(estuary_bottom, shelf_bottom; width, minimum_depth)
restrict_bottom!(shelf_bottom, estuary_bottom; width, minimum_depth)

levels = (shelf = (shelf_grid, shelf_bottom), estuary = (estuary_grid, estuary_bottom))

#####
##### Figure
#####

figure = Figure(size=(1400, 560))

for (n, name) in enumerate(keys(levels))
    grid, bottom = levels[name]
    λ, φ = λnodes(grid, Center()), φnodes(grid, Center())
    depth = -Array(interior(bottom, :, :, 1))
    depth[depth .≤ 0] .= NaN

    aspect = AxisAspect((λ[end] - λ[1]) * cosd((φ[1] + φ[end]) / 2) / (φ[end] - φ[1]))
    title = @sprintf("%s, %d × %d × %d", name, size(grid)...)
    axis = Axis(figure[1, n]; title, aspect, xlabel="longitude (°E)", ylabel = n == 1 ? "latitude (°N)" : "")
    depth_map = heatmap!(axis, λ, φ, depth; colormap=:deep, nan_color=:gray80)
    Colorbar(figure[2, n], depth_map; vertical=false, flipaxis=false, label="depth (m)")

    if n == 1
        (west, east), (south, north) = x_domain(estuary_grid), y_domain(estuary_grid)
        lines!(axis, [west, east, east, west, west], [south, south, north, north, south]; color=:red)
    end
end

outpath = joinpath(@__DIR__, "severn_bathymetry.png")
save(outpath, figure)
println("Saved to $outpath")

#####
##### Column volumes
#####

function column_volumes(grid, immersed_boundary)
    grid = ImmersedBoundaryGrid(grid, immersed_boundary)
    Nx, Ny, _ = size(grid)
    return [Azᶜᶜᶜ(i, j, 1, grid) * static_column_depthᶜᶜᵃ(i, j, grid) for i in 1:Nx, j in 1:Ny]
end

for (name, immersed_boundary) in (("partial cells", PartialCellBottom),
                                   ("partial cells of any thickness", bottom -> PartialCellBottom(bottom; minimum_fractional_cell_height=0)))
    parent_grid, parent_bottom = levels.shelf
    grid, bottom = levels.estuary
    i, j, _ = parent_indices(parent_grid, grid)
    inner = (1+width:length(i)-width, 1+width:length(j)-width)

    nest_volume = reshape(column_volumes(grid, immersed_boundary(bottom)),
                          refinement, length(i), refinement, length(j))
    nest_volume = dropdims(sum(nest_volume, dims=(1, 3)), dims=(1, 3))[inner...]
    parent_volume = column_volumes(parent_grid, immersed_boundary(parent_bottom))[i, j][inner...]

    mismatch = abs.(parent_volume .- nest_volume) ./ max.(parent_volume, nest_volume, 1)
    differing = mismatch .> 1e-6

    @printf("shelf ← estuary with %s: %d parent columns, %d wet\n",
            name, length(mismatch), count(>(0), parent_volume))
    @printf("    total volume: parent %.6e m³, nest %.6e m³, relative difference %.2e\n",
            sum(parent_volume), sum(nest_volume), (sum(parent_volume) - sum(nest_volume)) / sum(nest_volume))
    @printf("    columns differing by more than 1e-6: %d\n", count(differing))
end
