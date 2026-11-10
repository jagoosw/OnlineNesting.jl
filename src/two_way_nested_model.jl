struct TwoWayNestedModel{P, C, X, I, M, F, K, TS, A} <: AbstractModel{TS, A}
    parent :: P
    child :: C
    exchanger :: X
    refinement :: Int
    width :: Int
    feedback_indices :: I
    flux_mismatch :: M
    flux_callbacks :: F
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

# TODO: `time_step!(::EarthSystemModel, Δt; callbacks)` in NumericalEarth should pass `callbacks` to its ocean,
# as a model does to its stages; a level would then need no callback of its own.
add_ocean_callback!(level, callback) = nothing
add_ocean_callback!(level::EarthSystemModel, callback) = add_callback!(level.ocean, callback)
add_ocean_callback!(nest::TwoWayNestedModel, callback) = add_ocean_callback!(nest.parent, callback)

"""
$(TYPEDSIGNATURES)

Return a `child` nested in `parent`, with boundary conditions that read `exchanger`. The parent
is an ocean model or an `EarthSystemModel`; the child is either of these or another
`TwoWayNestedModel`. A time step advances the parent, then the child to the same time in as
many steps as there are child cells across a parent cell, then feeds the child's state back to
the parent in the parent cells the child covers, apart from its outermost `width` parent cells:
tracers and sea level become the volume- and area-weighted means over the child's wet cells,
and `u` and `v` become the child's volume transport through each parent face over the area of
the child's faces there. The parent cells around that region then take the difference between
the volume and tracer that the parent and the child carried through its edge during the step.
The parent outside the region and the child inside it then conserve volume together and, on a
z-star grid, tracer content, when `width` is at least 1 and nothing forces the child within the
region. The levels use a split Runge-Kutta time stepper. `Simulation` callbacks that act within
a time step are passed to the parent.
"""
function TwoWayNestedModel(parent::AbstractModel{TS, A}, child, exchanger; width) where {TS, A}
    grid = ocean_model(child).grid
    i, j, k = parent_indices(ocean_model(parent).grid, grid)
    refinement = size(grid, 1) ÷ length(i)
    feedback_indices = (i[1+width:end-width], j[1+width:end-width], k)

    parent_ocean = ocean_model(parent)
    flux_mismatch = (volume = CenterField(parent_ocean.grid), tracers = map(c -> CenterField(parent_ocean.grid), parent_ocean.tracers))
    first_cell = (first(i), first(j), first(k))

    parent_fluxes = (; flux_mismatch, feedback_indices, sign = 1, refinement = 1, first_cell = (1, 1, 1))
    child_fluxes = (; flux_mismatch, feedback_indices, sign = -1, refinement, first_cell)

    flux_callbacks = map((parent = parent_fluxes, child = child_fluxes)) do parameters
        Callback(accumulate_fluxes!; callsite = TracerTendencyCallsite(), parameters)
    end

    add_ocean_callback!(parent, flux_callbacks.parent)
    add_ocean_callback!(child, flux_callbacks.child)

    P, C, X, I, M, F, K = typeof.((parent, child, exchanger, feedback_indices, flux_mismatch, flux_callbacks, parent.clock))
    return TwoWayNestedModel{P, C, X, I, M, F, K, TS, A}(parent, child, exchanger, refinement, width, feedback_indices,
                                                         flux_mismatch, flux_callbacks, parent.clock)
end

function TimeSteppers.time_step!(nest::TwoWayNestedModel, Δt; callbacks=[])
    nest.clock.iteration == 0 && update_state!(nest, callbacks)

    time = nest.clock.time
    store_state!(nest.exchanger, 1)
    time_step!(nest.parent, Δt; callbacks = [callbacks..., nest.flux_callbacks.parent])
    store_state!(nest.exchanger, 2)

    for series in nest.exchanger.state
        series.times = range(time, nest.clock.time; length=2)
    end

    child_callbacks = [nest.flux_callbacks.child]

    for _ in 2:nest.refinement
        time_step!(nest.child, Δt / nest.refinement; callbacks = child_callbacks)
    end

    # the last step is the time remaining, so the child arrives at the parent's time exactly
    time_step!(nest.child, nest.clock.time - nest.child.clock.time; callbacks = child_callbacks)
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

    correct_fluxes!(parent, nest.flux_mismatch, i, j)

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

# Add to the parent cells around the feedback region, in the last stage of a time step of `model`, the volume
# and tracer that the step carries through the region's edge: out of those cells for the parent (`sign = 1`),
# into them for the child (`sign = -1`).
function accumulate_fluxes!(model, parameters)
    model.clock.stage == model.timestepper.Nstages || return nothing

    i, j, k = parameters.feedback_indices
    refinement = parameters.refinement
    Δt = convert(eltype(model.grid), parameters.sign * model.clock.last_stage_Δt)

    # the cells west and east, then south and north, with the shift from each to its face of the region
    accumulate_side!(model, parameters, (first(i)-1:first(i)-1, j, k), (1, 0, 0),  Δt, (1, refinement), tracer_flux_x)
    accumulate_side!(model, parameters, (last(i)+1:last(i)+1, j, k),   (0, 0, 0), -Δt, (1, refinement), tracer_flux_x)
    accumulate_side!(model, parameters, (i, first(j)-1:first(j)-1, k), (0, 1, 0),  Δt, (refinement, 1), tracer_flux_y)
    accumulate_side!(model, parameters, (i, last(j)+1:last(j)+1, k),   (0, 0, 0), -Δt, (refinement, 1), tracer_flux_y)

    return nothing
end

function accumulate_side!(model, parameters, cells, shift, Δt, block, flux)
    grid = model.grid
    arch = architecture(grid)
    parent_grid = parameters.flux_mismatch.volume.grid
    velocities = model.transport_velocities
    kernel_parameters = KernelParameters(cells...)
    first_face = parameters.first_cell .- shift
    state = (model.clock, fields(model), model.buoyancy)

    # the centered flux of a tracer equal to one is the volume flux
    launch!(arch, parent_grid, kernel_parameters, _accumulate_flux!, parameters.flux_mismatch.volume, grid, flux, Δt, block,
            parameters.refinement, first_face, Centered(eltype(grid)), velocities, OneField(), nothing, nothing, Val(0), state...)

    for (name, mismatch) in pairs(parameters.flux_mismatch.tracers)
        tracer_index = Val(findfirst(==(name), keys(model.tracers)))
        advecting_velocities = tracer_advecting_velocities(velocities, model.biogeochemistry, model.closure, model.closure_fields,
                                                           model.forcing[name], Val(name))

        launch!(arch, parent_grid, kernel_parameters, _accumulate_flux!, mismatch, grid, flux, Δt, block,
                parameters.refinement, first_face, model.advection[name], advecting_velocities, model.tracers[name],
                model.closure, model.closure_fields, tracer_index, state...)
    end

    return nothing
end

# Add to `mismatch` in the parent cell (i, j, k) `Δt` times the sum of `flux` over a `block` of faces of `grid`;
# `first_face` is the parent cell whose block starts at the first face of `grid`
@kernel function _accumulate_flux!(mismatch, grid, flux, Δt, block, refinement, first_face, args...)
    i, j, k = @index(Global, NTuple)
    i₁, j₁, k₁ = first_face
    k′ = k - k₁ + 1
    total_flux = zero(grid)

    for n in 1:block[2], m in 1:block[1]
        i′ = refinement * (i - i₁) + m
        j′ = refinement * (j - j₁) + n
        total_flux += flux(i′, j′, k′, grid, args...)
    end

    @inbounds mismatch[i, j, k] += Δt * total_flux
end

function correct_fluxes!(parent, flux_mismatch, i, j)
    grid = parent.grid
    arch = architecture(grid)
    columns = KernelParameters(first(i)-1:last(i)+1, first(j)-1:last(j)+1)

    for (name, mismatch) in pairs(flux_mismatch.tracers)
        launch!(arch, grid, columns, _correct_tracer!, parent.tracers[name], grid, mismatch, flux_mismatch.volume)
        fill!(mismatch, 0)
    end

    launch!(arch, grid, columns, _correct_sea_level!, parent.free_surface.displacement, grid, flux_mismatch.volume)
    fill!(flux_mismatch.volume, 0)

    return nothing
end

# The tracer of a column whose layers take the volumes `volume_mismatch` and the tracer contents
# `tracer_mismatch` while its sea level rises by their total volume: a wet layer thickens in proportion to its
# volume and passes the rest of what it takes upward with the mean of its tracer and the next layer's; a land
# layer passes all of it.
@kernel function _correct_tracer!(tracer, grid, tracer_mismatch, volume_mismatch)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)
    column_volume = zero(grid)
    column_mismatch = zero(grid)

    for k in 1:Nz
        column_volume += wet_volume(i, j, k, grid)
        column_mismatch += @inbounds volume_mismatch[i, j, k]
    end

    stretching = ifelse(column_volume > 0, column_mismatch / column_volume, zero(grid))
    upward_transport = zero(grid)
    upward_content = zero(grid)

    for k in 1:Nz
        volume = wet_volume(i, j, k, grid)
        layer_tracer = @inbounds tracer[i, j, k]
        tracer_above = @inbounds tracer[i, j, min(k + 1, Nz)]
        volume_taken = upward_transport + @inbounds volume_mismatch[i, j, k]
        content_taken = upward_content + @inbounds tracer_mismatch[i, j, k]

        transport_above = volume_taken - volume * stretching
        content_above = ifelse(volume > 0, transport_above * (layer_tracer + tracer_above) / 2, content_taken)
        content_change = content_taken - content_above - layer_tracer * (volume_taken - transport_above)

        @inbounds tracer[i, j, k] = ifelse(volume > 0, layer_tracer + content_change / (volume * (1 + stretching)), layer_tracer)
        upward_transport = transport_above
        upward_content = content_above
    end
end

@kernel function _correct_sea_level!(η, grid, volume_mismatch)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)
    area = wet_surface_area(i, j, Nz + 1, grid)
    column_mismatch = zero(grid)

    for k in 1:Nz
        column_mismatch += @inbounds volume_mismatch[i, j, k]
    end

    @inbounds η[i, j, Nz+1] += ifelse(area > 0, column_mismatch / area, zero(grid))
end
