# EMODnet bathymetry on the shelf, approaches, bay and headland grids, reconciled
# between levels. Draws each level's depth with the outline of its nest, reports how well the column volumes of each parent and nest agree, and whether
# a stratified ocean stays at rest over each level's bottom, each level run alone with closed sides.

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: rnodes, static_column_depthᶜᶜᵃ, x_domain, y_domain, λnodes, φnodes
using Oceananigans.Operators: Azᶜᶜᶜ
using NumericalEarth
using NumericalEarth.DataWrangling: BoundingBox, Metadatum
using OnlineNesting
using CairoMakie
using Printf

refinement = 5
width = 3 # parent cells
minimum_depth = 10
z = vcat(-reverse(cumsum(0.5 .* 1.076 .^ (0:49))), 0) # 50 layers, 0.5 m thick at the surface, down to 250 m

shelf_grid = LatitudeLongitudeGrid(size=(220, 245, 50), longitude=(-10.5, -2.5), latitude=(49, 54.5), halo=(7, 7, 4); z)
approaches_grid = nested_grid(shelf_grid; longitude=(-6.4, -4.6), latitude=(51.3, 52.4), depth=150, refinement)
bay_grid = nested_grid(approaches_grid; longitude=(-5.50, -5.06), latitude=(51.68, 51.94), depth=95, refinement)
headland_grid = nested_grid(bay_grid; longitude=(-5.342, -5.244), latitude=(51.850, 51.912), depth=72, refinement)

function emodnet_bottom_height(grid; interpolation_passes, major_basins=Inf)
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(grid), dir=joinpath(dirname(@__DIR__), "data"))
    return regrid_bathymetry(grid, metadatum; minimum_depth, interpolation_passes, major_basins)
end

shelf_bottom = emodnet_bottom_height(shelf_grid; interpolation_passes=4, major_basins=1)
approaches_bottom = emodnet_bottom_height(approaches_grid; interpolation_passes=2)
bay_bottom = emodnet_bottom_height(bay_grid; interpolation_passes=1)
headland_bottom = emodnet_bottom_height(headland_grid; interpolation_passes=1)

blend_parent_bottom!(approaches_bottom, shelf_bottom; width, minimum_depth)
blend_parent_bottom!(bay_bottom, approaches_bottom; width, minimum_depth)
blend_parent_bottom!(headland_bottom, bay_bottom; width, minimum_depth)
restrict_bottom!(bay_bottom, headland_bottom; width, minimum_depth)
restrict_bottom!(approaches_bottom, bay_bottom; width, minimum_depth)
restrict_bottom!(shelf_bottom, approaches_bottom; width, minimum_depth)

levels = (shelf = (shelf_grid, shelf_bottom), approaches = (approaches_grid, approaches_bottom),
          bay = (bay_grid, bay_bottom), headland = (headland_grid, headland_bottom))

#####
##### Figure
#####

figure = Figure(size=(1700, 560))

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

    if n < length(levels)
        (west, east), (south, north) = x_domain(levels[n+1][1]), y_domain(levels[n+1][1])
        lines!(axis, [west, east, east, west, west], [south, south, north, north, south]; color=:red)
    end
end

save(joinpath(@__DIR__, "pembrokeshire_bathymetry.png"), figure)

#####
##### Column volumes
#####

function column_volumes(grid, immersed_boundary)
    grid = ImmersedBoundaryGrid(grid, immersed_boundary)
    Nx, Ny, _ = size(grid)
    return [Azᶜᶜᶜ(i, j, 1, grid) * static_column_depthᶜᶜᵃ(i, j, grid) for i in 1:Nx, j in 1:Ny]
end

immersed_boundaries = (("partial cells", PartialCellBottom),
                       ("partial cells of any thickness", bottom -> PartialCellBottom(bottom; minimum_fractional_cell_height=0)),
                       ("full cells", GridFittedBottom))

for (name, immersed_boundary) in immersed_boundaries, (parent, nest) in ((:shelf, :approaches), (:approaches, :bay), (:bay, :headland))
    (parent_grid, parent_bottom), (grid, bottom) = levels[parent], levels[nest]
    i, j, _ = parent_indices(parent_grid, grid)
    inner = (1+width:length(i)-width, 1+width:length(j)-width)

    nest_volume = reshape(column_volumes(grid, immersed_boundary(bottom)), refinement, length(i), refinement, length(j))
    nest_volume = dropdims(sum(nest_volume, dims=(1, 3)), dims=(1, 3))[inner...]
    parent_volume = column_volumes(parent_grid, immersed_boundary(parent_bottom))[i, j][inner...]
    parent_depth = -min.(0, Array(interior(parent_bottom, i, j, 1)))[inner...]

    mismatch = abs.(parent_volume .- nest_volume) ./ max.(parent_volume, nest_volume, 1)
    differing = mismatch .> 1e-6
    rounded = parent_depth .≤ minimum_depth

    @printf("%s ← %s with %s: %d parent columns in the feedback region, %d wet\n", parent, nest, name, length(mismatch), count(>(0), parent_volume))
    @printf("    total volume: parent %.6e m³, nest %.6e m³, relative difference %.2e\n",
            sum(parent_volume), sum(nest_volume), (sum(parent_volume) - sum(nest_volume)) / sum(nest_volume))
    @printf("    columns differing by more than 1e-6: %d, of which %d at or below the minimum depth\n",
            count(differing), count(differing .& rounded))
    @printf("    mean relative difference over the wet columns: %.2e\n", sum(mismatch) / count(>(0), max.(parent_volume, nest_volume)))
    @printf("    largest relative difference in the other columns: %.2e\n", maximum(mismatch[.!rounded]))
end

#####
##### Rest
#####

Δt = 2minutes
stop_time = 20minutes

for (n, name) in enumerate(keys(levels))
    grid, bottom_height = levels[name]
    closed_grid = LatitudeLongitudeGrid(size=size(grid), longitude=x_domain(grid), latitude=y_domain(grid),
                                        z=rnodes(grid, Face()), halo=(7, 7, 4))
    closed_grid = ImmersedBoundaryGrid(closed_grid, PartialCellBottom(Array(interior(bottom_height, :, :, 1))))

    model = HydrostaticFreeSurfaceModel(closed_grid;
                                        tracers = :T,
                                        buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(), constant_salinity=35),
                                        coriolis = HydrostaticSphericalCoriolis(),
                                        free_surface = SplitExplicitFreeSurface(closed_grid; substeps=40),
                                        timestepper = :SplitRungeKutta3)

    set!(model, T = (λ, φ, z) -> 10 + 8 * exp(z / 50))

    simulation = Simulation(model; Δt = Δt / refinement^(n - 1), stop_time)
    run!(simulation)

    u, v, w = model.velocities
    @printf("%s at rest after %s (%d steps): max|u| = %.2e, max|v| = %.2e, max|w| = %.2e m s⁻¹, max|η| = %.2e m\n",
            name, prettytime(stop_time), iteration(simulation), maximum(abs, u), maximum(abs, v), maximum(abs, w),
            maximum(abs, model.free_surface.displacement))
end
