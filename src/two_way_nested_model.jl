struct TwoWayNestedModel{P, C, X, K, TS, A} <: AbstractModel{TS, A}
    parent :: P
    child :: C
    exchanger :: X
    refinement :: Int
    clock :: K
end

"""
$(TYPEDSIGNATURES)

Return a `child` ocean model nested in `parent`, with boundary conditions that read
`exchanger`. A time step advances the parent, then the child in `refinement` steps to the same
time. `refinement` defaults to the ratio of the grid spacings.
"""
function TwoWayNestedModel(parent::AbstractModel{TS, A}, child, exchanger;
                           refinement = size(child.grid, 1) ÷
                                        length(first(parent_indices(parent.grid, child.grid)))) where {TS, A}

    P, C, X, K = typeof(parent), typeof(child), typeof(exchanger), typeof(parent.clock)
    return TwoWayNestedModel{P, C, X, K, TS, A}(parent, child, exchanger, refinement, parent.clock)
end

function TimeSteppers.time_step!(nest::TwoWayNestedModel, Δt; callbacks=[])
    time = nest.clock.time
    store_state!(nest.exchanger, 1)
    time_step!(nest.parent, Δt)
    store_state!(nest.exchanger, 2)

    for series in nest.exchanger.state
        series.times = range(time, nest.clock.time; length=2)
    end

    for _ in 1:nest.refinement
        time_step!(nest.child, Δt / nest.refinement)
    end

    return nothing
end

function TimeSteppers.update_state!(nest::TwoWayNestedModel, callbacks=[])
    update_state!(nest.parent)
    update_state!(nest.child)
    return nothing
end
