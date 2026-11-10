# EMODnet bathymetry on the shelf, approaches and bay grids, reconciled
# between levels. Reports how well the column volumes of each parent and nest agree, and whether
# a stratified ocean stays at rest over each level's bottom, each level run alone with closed sides.

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: rnodes, static_column_depthᶜᶜᵃ, x_domain, y_domain
using Oceananigans.Operators: Azᶜᶜᶜ
using NumericalEarth
using NumericalEarth.DataWrangling: BoundingBox, Metadatum
using OnlineNesting
using Printf

refinement = 5
width = 3 # parent cells
minimum_depth = 10
z = [-250 * (1 - k / 14)^1.7 for k in 0:14]

shelf_grid = LatitudeLongitudeGrid(size=(220, 245, 14), longitude=(-10.5, -2.5), latitude=(49, 54.5), halo=(7, 7, 4); z)
approaches_grid = nested_grid(shelf_grid; longitude=(-6.4, -4.6), latitude=(51.3, 52.4), depth=150, refinement)
bay_grid = nested_grid(approaches_grid; longitude=(-5.50, -5.06), latitude=(51.68, 51.94), depth=95, refinement)

function emodnet_bottom_height(grid; interpolation_passes, major_basins=Inf)
    metadatum = Metadatum(:bottom_height; dataset=EMODnet2024(), region=BoundingBox(grid))
    return regrid_bathymetry(grid, metadatum; minimum_depth, interpolation_passes, major_basins)
end

shelf_bottom = emodnet_bottom_height(shelf_grid; interpolation_passes=4, major_basins=1)
approaches_bottom = emodnet_bottom_height(approaches_grid; interpolation_passes=2)
bay_bottom = emodnet_bottom_height(bay_grid; interpolation_passes=1)

blend_parent_bottom!(approaches_bottom, shelf_bottom; width, minimum_depth)
blend_parent_bottom!(bay_bottom, approaches_bottom; width, minimum_depth)
restrict_bottom!(approaches_bottom, bay_bottom; width, minimum_depth)
restrict_bottom!(shelf_bottom, approaches_bottom; width, minimum_depth)

levels = (shelf = (shelf_grid, shelf_bottom), approaches = (approaches_grid, approaches_bottom), bay = (bay_grid, bay_bottom))

#####
##### Column volumes
#####

function column_volumes(grid, bottom_height, minimum_fractional_cell_height)
    grid = ImmersedBoundaryGrid(grid, PartialCellBottom(bottom_height; minimum_fractional_cell_height))
    Nx, Ny, _ = size(grid)
    return [Azᶜᶜᶜ(i, j, 1, grid) * static_column_depthᶜᶜᵃ(i, j, grid) for i in 1:Nx, j in 1:Ny]
end

for minimum_fractional_cell_height in (0.2, 0), (parent, nest) in ((:shelf, :approaches), (:approaches, :bay))
    (parent_grid, parent_bottom), (grid, bottom) = levels[parent], levels[nest]
    i, j, _ = parent_indices(parent_grid, grid)
    inner = (1+width:length(i)-width, 1+width:length(j)-width)

    nest_volume = reshape(column_volumes(grid, bottom, minimum_fractional_cell_height), refinement, length(i), refinement, length(j))
    nest_volume = dropdims(sum(nest_volume, dims=(1, 3)), dims=(1, 3))[inner...]
    parent_volume = column_volumes(parent_grid, parent_bottom, minimum_fractional_cell_height)[i, j][inner...]
    parent_depth = -min.(0, Array(interior(parent_bottom, i, j, 1)))[inner...]

    mismatch = abs.(parent_volume .- nest_volume) ./ max.(parent_volume, nest_volume, 1)
    differing = mismatch .> 1e-6
    rounded = parent_depth .≤ minimum_depth

    @printf("%s ← %s with a minimum fractional cell height of %.1f: %d parent columns in the feedback region, %d wet\n",
            parent, nest, minimum_fractional_cell_height, length(mismatch), count(>(0), parent_volume))
    @printf("    total volume: parent %.6e m³, nest %.6e m³, relative difference %.2e\n",
            sum(parent_volume), sum(nest_volume), (sum(parent_volume) - sum(nest_volume)) / sum(nest_volume))
    @printf("    columns differing by more than 1e-6: %d, of which %d at or below the minimum depth\n",
            count(differing), count(differing .& rounded))
    @printf("    largest relative difference in the other columns: %.2e\n", maximum(mismatch[.!rounded]))
end

#####
##### Rest
#####

Δt = 2minutes
stop_time = 1hour

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
