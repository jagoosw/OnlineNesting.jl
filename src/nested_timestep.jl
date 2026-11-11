"""
$(TYPEDSIGNATURES)

Compute the minimum horizontal grid spacing [m] for a `LatitudeLongitudeGrid`, accounting for
the convergence of meridians. The spacing is evaluated at the latitude farthest from the equator,
where cells are narrowest.
"""
function minimum_horizontal_spacing(grid::LatitudeLongitudeGrid)
    longitude_west, longitude_east = x_domain(grid)
    latitude_south, latitude_north = y_domain(grid)
    number_of_x_cells, number_of_y_cells, _ = size(grid)

    degrees_to_radians = π / 180
    zonal_cell_degrees = (longitude_east - longitude_west) / number_of_x_cells
    meridional_cell_degrees = (latitude_north - latitude_south) / number_of_y_cells

    extreme_latitude = max(abs(latitude_south), abs(latitude_north))
    zonal_spacing = grid.radius * zonal_cell_degrees * degrees_to_radians * cos(extreme_latitude * degrees_to_radians)
    meridional_spacing = grid.radius * meridional_cell_degrees * degrees_to_radians

    return min(zonal_spacing, meridional_spacing)
end

"""
$(TYPEDSIGNATURES)

Compute the minimum horizontal grid spacing [m] for a `RectilinearGrid`.
"""
function minimum_horizontal_spacing(grid::RectilinearGrid)
    x_west, x_east = x_domain(grid)
    y_south, y_north = y_domain(grid)
    number_of_x_cells, number_of_y_cells, _ = size(grid)

    return min((x_east - x_west) / number_of_x_cells,
               (y_north - y_south) / number_of_y_cells)
end

"""
$(TYPEDSIGNATURES)

Compute timesteps for a hierarchy of nested grids from the advective CFL condition. The `grids`
are ordered from coarsest (level 1) to finest, with a matching collection of `maximum_velocities`
[m/s] at each level.

The nesting constraint requires each child to take `refinement` steps per parent step, where the
refinement ratio is inferred from adjacent grids. The binding CFL constraint at any level
determines all timesteps.

Returns a `NamedTuple` with `timesteps` (a `Tuple`, one per level) and `binding_level`.

Keyword arguments:
- `cfl_number`: CFL number (default `1.7`, appropriate for RK3).

```jldoctest
using Oceananigans
using OnlineNesting

shelf = LatitudeLongitudeGrid(size=(8, 6, 4), longitude=(-10, -2), latitude=(49, 55), z=(-200, 0))
child = nested_grid(shelf; longitude=(-6.2, -3.9), latitude=(51, 53), refinement=5, depth=100)
result = redirect_stdout(devnull) do
    nested_timestep((shelf, child), (2.0, 1.5))
end
length(result.timesteps)

# output
2
```
"""
function nested_timestep(grids, maximum_velocities; cfl_number=1.7)
    number_of_levels = length(grids)

    horizontal_spacings = [minimum_horizontal_spacing(grids[i]) for i in 1:number_of_levels]
    cfl_limited_timesteps = [cfl_number * horizontal_spacings[i] / maximum_velocities[i]
                             for i in 1:number_of_levels]

    refinement_ratios = Vector{Int}(undef, number_of_levels - 1)
    for level in 2:number_of_levels
        parent_cell_extent = (x_domain(grids[level - 1])[2] - x_domain(grids[level - 1])[1]) / size(grids[level - 1], 1)
        child_cell_extent  = (x_domain(grids[level])[2] - x_domain(grids[level])[1]) / size(grids[level], 1)
        refinement_ratios[level - 1] = round(Int, parent_cell_extent / child_cell_extent)
    end

    implied_coarsest_timesteps = [cfl_limited_timesteps[i] * prod(refinement_ratios[1:i-1])
                                  for i in 1:number_of_levels]

    binding_level = argmin(implied_coarsest_timesteps)
    coarsest_timestep = implied_coarsest_timesteps[binding_level]

    timesteps = Tuple(coarsest_timestep / prod(refinement_ratios[1:i-1]) for i in 1:number_of_levels)

    println("Nested timestep summary ($number_of_levels levels, CFL = $cfl_number):")
    println("─" ^ 78)
    println("  ", rpad("Level", 8), "  ", lpad("Nx × Ny", 10), "  ",
            lpad("Δx [m]", 10), "  ", lpad("U_max", 10), "  ",
            lpad("CFL Δt [s]", 12), "  ", lpad("Δt [s]", 10))
    println("─" ^ 78)
    for i in 1:number_of_levels
        number_of_x_cells, number_of_y_cells, _ = size(grids[i])
        grid_size_string = "$(number_of_x_cells)×$(number_of_y_cells)"
        marker = i == binding_level ? " ←" : ""
        println("  ", rpad(i, 8), "  ", lpad(grid_size_string, 10), "  ",
                lpad(round(Int, horizontal_spacings[i]), 10), "  ",
                lpad(round(maximum_velocities[i], digits=2), 10), "  ",
                lpad(round(cfl_limited_timesteps[i], digits=1), 12), "  ",
                lpad(round(timesteps[i], digits=1), 10), marker)
    end
    println("─" ^ 78)
    println("Binding constraint: level $binding_level")

    return (; timesteps, binding_level)
end
