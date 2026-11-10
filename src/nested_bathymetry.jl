sea_floor_height(bottom_height) = clamp.(Array(interior(bottom_height, :, :, 1)), -bottom_height.grid.Lz, 0)

# Without wetting and drying a column is land or at least `minimum_depth` deep, whichever is nearer
function limit_depth(height, minimum_depth)
    limited_height = min(height, -minimum_depth)
    return ifelse(height > -minimum_depth / 2, zero(limited_height), limited_height)
end

"""
$(TYPEDSIGNATURES)

Refer `bottom_height` [m], measured from the lowest astronomical tide, to mean sea level by
adding `lowest_tide`, the height [m] of the lowest astronomical tide relative to mean sea level
(negative): a number, a function of longitude and latitude, or a field, which is interpolated to
the grid of `bottom_height`. Land is set to zero height; columns left shallower than half of
`minimum_depth` [m] become land and the other wet columns are at least `minimum_depth` deep.
"""
function refer_to_mean_sea_level!(bottom_height, lowest_tide; minimum_depth = 0)
    datum_height = Field{Center, Center, Nothing}(bottom_height.grid)
    set!(datum_height, lowest_tide)

    height = Array(interior(bottom_height, :, :, 1)) .+ Array(interior(datum_height, :, :, 1))
    set!(bottom_height, limit_depth.(height, minimum_depth))
    fill_halo_regions!(bottom_height)

    return bottom_height
end

"""
$(TYPEDSIGNATURES)

Blend the `bottom_height` [m] of a nest toward `parent_bottom_height`, the bottom height on its
parent grid, over the outermost `width` parent cells of the nest. The nest takes its parent's
bottom along its boundary and keeps its own `width` parent cells inward, with a linear ramp
between. Land is set to zero height and the sea floor is kept within the depth of the nest grid.
Columns left shallower than half of `minimum_depth` [m] become land and the other wet columns
are at least `minimum_depth` deep.
"""
function blend_parent_bottom!(bottom_height, parent_bottom_height; width, minimum_depth = 0)
    grid = bottom_height.grid
    i, j, _ = parent_indices(parent_bottom_height.grid, grid)
    refinement = size(grid, 1) ÷ length(i)

    parent_height = Field{Center, Center, Nothing}(grid)
    set!(parent_height, repeat(sea_floor_height(parent_bottom_height)[i, j], inner=(refinement, refinement)))
    set!(bottom_height, sea_floor_height(bottom_height))
    blend_parent_terrain!(bottom_height, parent_height; width = refinement * width)

    set!(bottom_height, limit_depth.(sea_floor_height(bottom_height), minimum_depth))
    fill_halo_regions!(bottom_height)

    return bottom_height
end

"""
$(TYPEDSIGNATURES)

Replace `parent_bottom_height` [m] by the area-weighted mean of the nest's `bottom_height`, with
land at zero height and the sea floor kept within the depth of the nest grid, in the parent
cells the nest covers apart from its outermost `width` parent cells. A parent column then holds
the volume of the nest columns within it, unless it is shallower than `minimum_depth` [m]:
columns shallower than half of `minimum_depth` become land and the others are at least
`minimum_depth` deep.
"""
function restrict_bottom!(parent_bottom_height, bottom_height; width, minimum_depth = 0)
    grid = bottom_height.grid
    i, j, _ = parent_indices(parent_bottom_height.grid, grid)
    Nx, Ny = length(i), length(j)
    refinement = size(grid, 1) ÷ Nx

    λ, φ, _ = cpu_faces(grid)
    area = diff(λ) .* diff(sind.(φ))'
    block_sum(a) = dropdims(sum(reshape(a, refinement, Nx, refinement, Ny), dims=(1, 3)), dims=(1, 3))
    mean_height = block_sum(area .* sea_floor_height(bottom_height)) ./ block_sum(area)

    inner = (1+width:Nx-width, 1+width:Ny-width)
    parent_height = Array(interior(parent_bottom_height, :, :, 1))
    parent_height[i[inner[1]], j[inner[2]]] .= limit_depth.(mean_height[inner...], minimum_depth)
    set!(parent_bottom_height, parent_height)
    fill_halo_regions!(parent_bottom_height)

    return parent_bottom_height
end
