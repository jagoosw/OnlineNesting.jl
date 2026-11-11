# EMODnet bathymetry on three nested levels (NW shelf → Bristol Channel →
# Severn Estuary), reconciled between levels. Draws each level's depth with the
# outline of the next nest, and reports column volume agreement.

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

shelf_grid = LatitudeLongitudeGrid(size=(340, 320, Nz),
    longitude=(-12.0, 5.0), latitude=(46.0, 62.0), halo=(7, 7, 4); z=z_faces)

bristol_grid = nested_grid(shelf_grid;
    longitude=(-5.5, -2.0), latitude=(50.8, 52.2), depth=250, refinement)

severn_grid = nested_grid(bristol_grid;
    longitude=(-3.5, -2.5), latitude=(51.1, 51.75), depth=100, refinement)

data_dir = joinpath(dirname(@__DIR__), "data")
tpxo_dir = get(ENV, "TPXO_DIR", joinpath(homedir(), "Downloads", "TPXO10_atlas_v2_nc"))

tpxo = TPXO10Atlas(; dir=tpxo_dir)
lat_field = lowest_astronomical_tide(tpxo, shelf_grid)
λ_west = x_domain(lat_field.grid)[1]
lowest_tide(λ, φ) = interpolate((mod(λ - λ_west, 360) + λ_west, φ), lat_field)

function emodnet_bottom_height(grid; interpolation_passes, major_basins=Inf, lat_correction=lowest_tide)
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(grid), dir=data_dir)
    bottom_height = regrid_bathymetry(grid, metadatum; minimum_depth, interpolation_passes, major_basins)
    return refer_to_mean_sea_level!(bottom_height, lat_correction; minimum_depth)
end

println("Downloading and regridding bathymetry...")
shelf_bottom   = emodnet_bottom_height(shelf_grid;   interpolation_passes=4, major_basins=1, lat_correction=0)
bristol_bottom = emodnet_bottom_height(bristol_grid;  interpolation_passes=2)
severn_bottom  = emodnet_bottom_height(severn_grid;   interpolation_passes=1)

blend_parent_bottom!(severn_bottom,  bristol_bottom; width, minimum_depth)
restrict_bottom!(bristol_bottom, severn_bottom;  width, minimum_depth)
restrict_bottom!(shelf_bottom,   bristol_bottom; width, minimum_depth)

levels = (shelf=(shelf_grid, shelf_bottom), bristol=(bristol_grid, bristol_bottom), severn=(severn_grid, severn_bottom))
child_grids = (bristol_grid, severn_grid, nothing)
color_ranges = ((0, 300), (0, 150), (0, 80))

#####
##### Figure
#####

figure = Figure(size=(1800, 700))

for (n, name) in enumerate(keys(levels))
    grid, bottom = levels[name]
    Nx, Ny, Nz_grid = size(grid)
    Δλ = (x_domain(grid)[2] - x_domain(grid)[1]) / Nx
    ϕ = sum(y_domain(grid)) / 2
    Δx = Δλ * π/180 * 6371e3 * cosd(ϕ)
    wet_fraction = count(Array(interior(bottom, :, :, 1)) .< 0) / (Nx * Ny)

    λ = Array(λnodes(grid, Center()))
    φ = Array(φnodes(grid, Center()))
    depth = -Array(interior(bottom, :, :, 1))
    depth[depth .≤ 0] .= NaN

    aspect = AxisAspect((λ[end] - λ[1]) * cosd((φ[1] + φ[end]) / 2) / (φ[end] - φ[1]))
    title = @sprintf("%s, %d × %d × %d (Δx ≈ %.0f m, %.0f%% wet)", name, Nx, Ny, Nz_grid, Δx, 100 * wet_fraction)
    axis = Axis(figure[1, n]; title, aspect, xlabel="longitude (°E)", ylabel = n == 1 ? "latitude (°N)" : "")
    depth_map = heatmap!(axis, λ, φ, depth; colormap=:deep, nan_color=:gray80, colorrange=color_ranges[n])
    Colorbar(figure[2, n], depth_map; vertical=false, flipaxis=false, label="depth (m)")

    child = child_grids[n]
    if !isnothing(child)
        (west, east), (south, north) = x_domain(child), y_domain(child)
        lines!(axis, [west, east, east, west, west], [south, south, north, north, south]; color=:red, linewidth=2)
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
    for (parent_name, child_name) in ((:shelf, :bristol), (:bristol, :severn))
        parent_grid, parent_bottom = levels[parent_name]
        grid, bottom = levels[child_name]
        i, j, _ = parent_indices(parent_grid, grid)
        inner = (1+width:length(i)-width, 1+width:length(j)-width)

        nest_volume = reshape(column_volumes(grid, immersed_boundary(bottom)),
                              refinement, length(i), refinement, length(j))
        nest_volume = dropdims(sum(nest_volume, dims=(1, 3)), dims=(1, 3))[inner...]
        parent_volume = column_volumes(parent_grid, immersed_boundary(parent_bottom))[i, j][inner...]

        mismatch = abs.(parent_volume .- nest_volume) ./ max.(parent_volume, nest_volume, 1)
        differing = mismatch .> 1e-6

        @printf("%s ← %s with %s: %d parent columns, %d wet\n",
                parent_name, child_name, name, length(mismatch), count(>(0), parent_volume))
        @printf("    total volume: parent %.6e m³, nest %.6e m³, relative difference %.2e\n",
                sum(parent_volume), sum(nest_volume), (sum(parent_volume) - sum(nest_volume)) / sum(nest_volume))
        @printf("    columns differing by more than 1e-6: %d\n", count(differing))
    end
end
