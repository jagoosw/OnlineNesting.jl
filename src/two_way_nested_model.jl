struct TwoWayNestedModel{P, C, X, I, K, TS, A} <: AbstractModel{TS, A}
    parent :: P
    child :: C
    exchanger :: X
    refinement :: Int
    width :: Int
    feedback_indices :: I
    clock :: K
end

ocean_model(model) = model
ocean_model(model::EarthSystemModel) = model.ocean.model
ocean_model(nest::TwoWayNestedModel) = ocean_model(nest.parent)

update_level_state!(level, callbacks=[]) = update_state!(level, callbacks)

function update_level_state!(level::EarthSystemModel, callbacks=[])
    update_state!(level.ocean.model, callbacks)
    update_state!(level)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return a `child` nested in `parent`, with boundary conditions that read `exchanger`. The parent
is an ocean model or an `EarthSystemModel`; the child is either of these or another
`TwoWayNestedModel`. A time step advances the parent, then the child to the same time in as
many steps as there are child cells across a parent cell, then feeds the child's state back to
the parent in the parent cells the child covers, apart from its outermost `width` parent cells:
tracers and sea level become the volume- and area-weighted means over the child's wet cells,
and `u` and `v` become the child's volume transport through each parent face over the area of
the child's faces there. `Simulation` callbacks that act within a time step are passed to the
parent.
"""
function TwoWayNestedModel(parent::AbstractModel{TS, A}, child, exchanger; width) where {TS, A}
    grid = ocean_model(child).grid
    i, j, k = parent_indices(ocean_model(parent).grid, grid)
    refinement = size(grid, 1) ÷ length(i)
    feedback_indices = (i[1+width:end-width], j[1+width:end-width], k)

    P, C, X, I, K = typeof(parent), typeof(child), typeof(exchanger), typeof(feedback_indices), typeof(parent.clock)
    return TwoWayNestedModel{P, C, X, I, K, TS, A}(parent, child, exchanger, refinement, width, feedback_indices, parent.clock)
end

function TimeSteppers.time_step!(nest::TwoWayNestedModel, Δt; callbacks=[])
    nest.clock.iteration == 0 && update_state!(nest, callbacks)

    time = nest.clock.time
    store_state!(nest.exchanger, 1)
    time_step!(nest.parent, Δt; callbacks)
    store_state!(nest.exchanger, 2)

    for series in nest.exchanger.state
        series.times = range(time, nest.clock.time; length=2)
    end

    for _ in 2:nest.refinement
        time_step!(nest.child, Δt / nest.refinement)
    end

    # the last step is the time remaining, so the child arrives at the parent's time exactly
    time_step!(nest.child, nest.clock.time - nest.child.clock.time)
    feed_back!(nest, callbacks)

    return nothing
end

function TimeSteppers.update_state!(nest::TwoWayNestedModel, callbacks=[])
    update_level_state!(nest.child)
    feed_back!(nest, callbacks)
    return nothing
end

function feed_back!(nest, callbacks)
    parent, child = ocean_model(nest.parent), ocean_model(nest.child)
    arch = architecture(parent.grid)
    refinement = nest.refinement
    i, j, k = nest.feedback_indices
    first_cell = (first(i) - nest.width, first(j) - nest.width, first(k))
    cells = (refinement, refinement)

    for (name, parent_tracer) in pairs(parent.tracers)
        launch!(arch, parent.grid, KernelParameters(i, j, k), _restrict!, parent_tracer, child.tracers[name],
                child.grid, wet_volume, cells, refinement, first_cell)
    end

    launch!(arch, parent.grid, KernelParameters(i, j, last(k)+1:last(k)+1), _restrict!, parent.free_surface.displacement,
            child.free_surface.displacement, child.grid, wet_surface_area, cells, refinement, first_cell)

    # the faces on both sides of every cell of the feedback region
    launch!(arch, parent.grid, KernelParameters(first(i):last(i)+1, j, k), _restrict!, parent.velocities.u,
            child.velocities.u, child.grid, Axᶠᶜᶜ, (1, refinement), refinement, first_cell)

    launch!(arch, parent.grid, KernelParameters(i, first(j):last(j)+1, k), _restrict!, parent.velocities.v,
            child.velocities.v, child.grid, Ayᶜᶠᶜ, (refinement, 1), refinement, first_cell)

    reconcile_state!(parent)
    update_level_state!(nest.parent, callbacks)

    return nothing
end

@inline wet_volume(i, j, k, grid) = Vᶜᶜᶜ(i, j, k, grid) * !inactive_cell(i, j, k, grid)
@inline wet_surface_area(i, j, k, grid) = Azᶜᶜᶠ(i, j, k, grid) * !inactive_cell(i, j, k - 1, grid)

# The mean of `field`, weighted by `weight`, over the `block` of nodes of `grid` on the parent node (i, j, k).
# `first_cell` is the first parent cell `grid` covers.
@kernel function _restrict!(parent_field, field, grid, weight, block, refinement, first_cell)
    i, j, k = @index(Global, NTuple)
    i₁, j₁, k₁ = first_cell
    k′ = k - k₁ + 1
    total_weight = zero(grid)
    content = zero(grid)

    for n in 1:block[2], m in 1:block[1]
        i′ = refinement * (i - i₁) + m
        j′ = refinement * (j - j₁) + n
        node_weight = weight(i′, j′, k′, grid)
        total_weight += node_weight
        content += @inbounds node_weight * field[i′, j′, k′]
    end

    @inbounds parent_field[i, j, k] = ifelse(total_weight > 0, content / total_weight, parent_field[i, j, k])
end
