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
many steps as there are child cells across a parent cell, then replaces the parent's tracers by
the volume-weighted mean of the child's in the parent cells the child covers, apart from its
outermost `width` parent cells. `Simulation` callbacks that act within a time step are passed
to the parent.
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
    i, j, k = nest.feedback_indices

    for (name, parent_tracer) in pairs(parent.tracers)
        launch!(architecture(parent.grid), parent.grid, KernelParameters(i, j, k), _restrict_tracer!, parent_tracer,
                child.tracers[name], child.grid, nest.refinement, first(i) - nest.width, first(j) - nest.width, first(k))
    end

    update_level_state!(nest.parent, callbacks)

    return nothing
end

# (i′, j′, k′) are the cells of `grid` within the parent cell (i, j, k); (i₁, j₁, k₁) is the first parent cell `grid` covers
@kernel function _restrict_tracer!(parent_tracer, tracer, grid, refinement, i₁, j₁, k₁)
    i, j, k = @index(Global, NTuple)
    k′ = k - k₁ + 1
    volume = zero(grid)
    content = zero(grid)

    for n in 1:refinement, m in 1:refinement
        i′ = refinement * (i - i₁) + m
        j′ = refinement * (j - j₁) + n
        V = Vᶜᶜᶜ(i′, j′, k′, grid) * !inactive_cell(i′, j′, k′, grid)
        volume += V
        content += @inbounds V * tracer[i′, j′, k′]
    end

    @inbounds parent_tracer[i, j, k] = ifelse(volume > 0, content / volume, parent_tracer[i, j, k])
end
