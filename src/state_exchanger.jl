struct StateExchanger{F, S, L}
    parent_fields :: F
    state :: S
    land_fill :: L
end

# The sea level is interpolated horizontally, at the location of a barotropic transport
series_location(field) = instantiated_location(field)
series_location(::Field{Center, Center, Face}) = (Center(), Center(), nothing)

take_unfilled!(donor, neighbor) = donor .= ifelse.(donor .== 0, neighbor, donor)

# Indices into `parent(field)` of the land nodes of `field` and, for each, of the nearest sea node on its level
function land_fill(field)
    grid = field.grid
    LX, LY, _ = location(field)
    land_operation = KernelFunctionOperation{LX, LY, Center}(inactive_node, grid, LX(), LY(), Center())
    land = Array(interior(compute!(Field(land_operation)))) .== 1
    land = size(field, 3) == 1 ? land[:, :, end:end] : land
    Nx, Ny, _ = size(land)

    donor = ifelse.(land, 0, LinearIndices(land))
    unfilled = typemax(Int)

    while 0 < count(iszero, donor) < unfilled
        unfilled = count(iszero, donor)
        take_unfilled!(view(donor, 2:Nx, :, :), view(donor, 1:Nx-1, :, :))
        take_unfilled!(view(donor, 1:Nx-1, :, :), view(donor, 2:Nx, :, :))
        take_unfilled!(view(donor, :, 2:Ny, :), view(donor, :, 1:Ny-1, :))
        take_unfilled!(view(donor, :, 1:Ny-1, :), view(donor, :, 2:Ny, :))
    end

    indices = LinearIndices(parent(field))[parentindices(interior(field))...]

    return (land = on_architecture(architecture(grid), indices[land]),
            donor = on_architecture(architecture(grid), indices[donor[land]]))
end

"""
$(TYPEDSIGNATURES)

Return an exchanger whose `state` holds a two-level `FieldTimeSeries` on the parent grid for
each prognostic field of the ocean model `parent`: `u`, `v`, `η`, the barotropic transports
`U` and `V`, and the tracers. The levels are the parent state at the start and at the end of
its latest time step, for a nest to interpolate between. Land nodes hold the value of the
nearest sea node on their level.
"""
function StateExchanger(parent)
    time = parent.clock.time
    parent_fields = prognostic_fields(parent)

    state = map(parent_fields) do field
        FieldTimeSeries(series_location(field), field.grid, range(time, time + 1; length=2))
    end

    return StateExchanger(parent_fields, state, map(land_fill, parent_fields))
end

function store_state!(exchanger, level)
    for (series, field, fill) in zip(exchanger.state, exchanger.parent_fields, exchanger.land_fill)
        data = view(parent(series), :, :, :, level)
        data .= parent(field)
        data[fill.land] = data[fill.donor]
    end

    return nothing
end

tracer_names(state) = filter(∉((:u, :v, :η, :U, :V)), keys(state))

"""
$(TYPEDSIGNATURES)

Return the lateral boundary conditions of an ocean model on the nest `grid` whose open sides
take the parent state held by `exchanger`:

- `U`, `V`: gravity-wave radiation toward the parent's barotropic transport and sea level;
- `u`, `v`: the parent's velocity, normal and tangential to each side;
- tracers: a `TracerReservoir` with the parent's tracer as its exterior value.

They are passed to the model as `boundary_conditions`.
"""
function nested_boundary_conditions(exchanger, grid)
    state = exchanger.state
    tracers = tracer_names(state)

    normal, tangential = NormalFlowBoundaryCondition, ValueBoundaryCondition
    reservoir(tracer) = ValueBoundaryCondition(tracer; scheme=TracerReservoir())

    bc_types = (u = (west=normal, east=normal, south=tangential, north=tangential),
                v = (west=tangential, east=tangential, south=normal, north=normal),
                NamedTuple{tracers}(map(Returns(reservoir), tracers))...)

    boundary_conditions = parent_boundary_conditions(grid; variables=state[(:u, :v, tracers...)], bc_types)

    # the elements of a tuple condition are tagged with their side here
    function radiation(transport, location, dimension, side)
        exterior = map((transport, state.η)) do series
            regularize_boundary_condition(Interpolated(series), grid, location, dimension, side)
        end
        return GravityWaveRadiationBoundaryCondition(exterior)
    end

    U = FieldBoundaryConditions(west = radiation(state.U, (Face(), Center(), nothing), 1, LeftBoundary),
                                east = radiation(state.U, (Face(), Center(), nothing), 1, RightBoundary))

    V = FieldBoundaryConditions(south = radiation(state.V, (Center(), Face(), nothing), 2, LeftBoundary),
                                north = radiation(state.V, (Center(), Face(), nothing), 2, RightBoundary))

    return merge(boundary_conditions, (; U, V))
end

"""
$(TYPEDSIGNATURES)

Return forcings that relax the velocities and tracers of an ocean model on the nest `grid`
toward the parent state held by `exchanger`, at `rate` [s⁻¹] on the open sides, decreasing
linearly to zero `width` cells into the nest.
"""
function nested_forcings(exchanger, grid; rate, width)
    Nx, Ny, _ = size(grid)
    λᶠ, φᶠ, _ = cpu_faces(grid)
    λ₁, λ₂ = extrema(λᶠ)
    φ₁, φ₂ = extrema(φᶠ)
    Lλ = width * (λ₂ - λ₁) / Nx
    Lφ = width * (φ₂ - φ₁) / Ny

    mask(λ, φ, z) = max(0, 1 - min((λ - λ₁) / Lλ, (λ₂ - λ) / Lλ, (φ - φ₁) / Lφ, (φ₂ - φ) / Lφ))
    state = exchanger.state

    return parent_forcings(; variables=state[(:u, :v, tracer_names(state)...)], rate, mask)
end
