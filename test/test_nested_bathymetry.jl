using Test
using Oceananigans
using Oceananigans.Grids: static_column_depthᶜᶜᵃ
using Oceananigans.Operators: Azᶜᶜᶜ
using OnlineNesting

parent_grid = LatitudeLongitudeGrid(size=(12, 10, 4), longitude=(-10, -2), latitude=(49, 55), z=[-100, -60, -30, -10, 0])
minimum_depth = 10
width = 1

# a sloping sea floor with land in the east; the nest resolves a bank and a shoaling, curved coast the parent does not
parent_bottom(λ, φ) = λ > -4.5 ? 20 : -40 + 5 * (λ + 6) - 3 * (φ - 52)
bank(λ, φ) = 25 * exp(-((λ + 6.5)^2 + (φ - 52)^2) / 0.1)
coast(φ) = -5.1 + 0.3 * sinpi(φ / 2)
nest_bottom(λ, φ) = λ > coast(φ) ? 35 : max(parent_bottom(λ, φ) + bank(λ, φ), min(-minimum_depth, 30 * (λ - coast(φ))))

function bottom_height_field(grid, bottom)
    bottom_height = Field{Center, Center, Nothing}(grid)
    set!(bottom_height, bottom)
    return bottom_height
end

heights(bottom_height) = Array(interior(bottom_height, :, :, 1))

function column_volumes(grid, bottom_height)
    grid = ImmersedBoundaryGrid(grid, PartialCellBottom(bottom_height; minimum_fractional_cell_height=0))
    Nx, Ny, _ = size(grid)
    return [Azᶜᶜᶜ(i, j, 1, grid) * static_column_depthᶜᶜᵃ(i, j, grid) for i in 1:Nx, j in 1:Ny]
end

@testset "nested bathymetry, refinement $refinement" for refinement in (2, 5)
    grid = nested_grid(parent_grid; longitude=(-8, -3.4), latitude=(50.2, 53.8), refinement)
    i, j, _ = parent_indices(parent_grid, grid)
    Nx, Ny, _ = size(grid)
    rim = refinement * width
    blocks(a) = dropdims(sum(reshape(a, refinement, length(i), refinement, length(j)), dims=(1, 3)), dims=(1, 3))

    parent_bottom_height = bottom_height_field(parent_grid, parent_bottom)
    bottom_height = bottom_height_field(grid, nest_bottom)
    original_parent_height = heights(parent_bottom_height)
    parent_height, height = min.(0, original_parent_height), min.(0, heights(bottom_height))

    blend_parent_bottom!(bottom_height, parent_bottom_height; width, minimum_depth)
    blended = heights(bottom_height)
    parent_on_nest = repeat(parent_height[i, j], inner=(refinement, refinement))
    boundary = [n ∈ (1, Nx) || m ∈ (1, Ny) for n in 1:Nx, m in 1:Ny]

    @test blended[boundary] == parent_on_nest[boundary]
    @test blended[1+rim:Nx-rim, 1+rim:Ny-rim] == height[1+rim:Nx-rim, 1+rim:Ny-rim]
    @test all(h -> h == 0 || h ≤ -minimum_depth, blended)

    # a linear ramp from the parent's bottom to the nest's, along a row in deep water
    m = Ny ÷ 2
    ramp = (0:rim) ./ rim
    @test blended[1:rim+1, m] ≈ ramp .* height[1:rim+1, m] .+ (1 .- ramp) .* parent_on_nest[1:rim+1, m]
    @test parent_on_nest[1, m] != height[1, m]

    restrict_bottom!(parent_bottom_height, bottom_height; width, minimum_depth)
    restricted = heights(parent_bottom_height)
    inner = (i[1+width:end-width], j[1+width:end-width])
    outside = trues(size(restricted))
    outside[inner...] .= false

    @test restricted[outside] == original_parent_height[outside]
    @test all(h -> h == 0 || h ≤ -minimum_depth, restricted[inner...])

    parent_volume = column_volumes(parent_grid, parent_bottom_height)[inner...]
    nest_volume = blocks(column_volumes(grid, bottom_height))[1+width:end-width, 1+width:end-width]
    deep = restricted[inner...] .< -minimum_depth
    @test count(deep) > length(deep) / 2
    @test parent_volume[deep] ≈ nest_volume[deep]

    # shallower columns are rounded to land or to the minimum depth, whichever is nearer
    for (nest_height, rounded_height) in ((-4, 0), (-6, -minimum_depth), (-30, -30))
        restrict_bottom!(parent_bottom_height, bottom_height_field(grid, nest_height); width, minimum_depth)
        @test all(heights(parent_bottom_height)[inner...] .≈ rounded_height)
    end
end

@testset "refer_to_mean_sea_level!" begin
    # heights from the lowest tide, across a beach: deep water, shallows and land
    beach(λ, φ) = 5 * (λ + 6) - 2
    lowest_tide(λ, φ) = -3 - (φ - 52) / 2
    tide_grid = LatitudeLongitudeGrid(size=(4, 5, 1), longitude=(-11, -1), latitude=(48, 56), z=(-1, 0))

    for tide in (lowest_tide, bottom_height_field(tide_grid, lowest_tide), -3)
        bottom_height = bottom_height_field(parent_grid, beach)
        refer_to_mean_sea_level!(bottom_height, tide; minimum_depth)

        expected = tide isa Number ? (λ, φ) -> beach(λ, φ) + tide : (λ, φ) -> beach(λ, φ) + lowest_tide(λ, φ)
        expected = heights(bottom_height_field(parent_grid, expected))
        referred = heights(bottom_height)

        deep, land = expected .< -minimum_depth, expected .> -minimum_depth / 2
        @test any(deep) && any(land) && any(.!deep .& .!land)
        @test referred[deep] ≈ expected[deep]
        @test all(referred[land] .== 0)
        @test all(referred[.!deep .& .!land] .== -minimum_depth)
    end
end
