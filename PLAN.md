# OnlineNesting.jl: plan for two-way, multi-level nesting

Target: a chain of ocean models that run together, coarse to fine, where each nest takes its
lateral boundary state from its parent and returns its full state (tracers, velocities,
barotropic transport and sea level) to the parent. The outermost level is any NumericalEarth
(NE) model: a global `EarthSystemModel`, or a regional one with open boundaries. The test
case is three levels: NW European shelf → Pembrokeshire approaches → St Brides Bay.

Versions this plan was checked against: Oceananigans 0.117.2, NumericalEarth 0.9.2, Julia 1.12.
Nothing in this plan has been run; every check below is still to be done.

## 1. Decisions

Settled:
- **Everything is two-way.** Tracers, `u`, `v`, the barotropic transport and `η` are all fed
  back. The feedback is built up field by field across milestones 4 and 5, but all of it is
  core scope.

Also settled: D1–D6 below, each confirmed at its default on 2026-11-10. The alternatives are
kept for the record.

| # | Decision | Chosen | Alternative not taken |
|---|----------|---------|-------------|
| D1 | Grid relationship | Child faces coincide with parent faces, integer refinement ratio, same grid type (`LatitudeLongitudeGrid`), shared vertical faces. Restriction is then an exact block sum. | General overlapping grids through ConservativeRegridding.jl (already an Oceananigans extension). More flexible, CPU-built weights, an extra dependency, and unverified for partial cells and face-located fields. |
| D2 | Conservation target | Restriction first, with the drift measured (milestones 4–5). Flux correction at the interface afterwards (milestone 7), which makes composite volume and tracer content conserved to round-off. | Stop after milestone 5 and accept the measured drift. |
| D3 | Surface forcing of the inner levels | Every level is its own NE model (an `EarthSystemModel`), so each computes its surface fluxes on its own grid. Each level gets its own prescribed-atmosphere instance over the same data, because `EarthSystemModel.time_step!` advances its atmosphere's clock on every call. | Inner levels are bare `HydrostaticFreeSurfaceModel`s that take fluxes interpolated from the parent. |
| D4 | Where the code lives | In this package, built on `NumericalEarth.NestedModels`. Pieces that prove general go upstream afterwards. | Develop directly on a NumericalEarth branch. |
| D5 | Which rule set governs | NumericalEarth's `AGENTS.md` and `.claude/rules`, including the restraint rules, since this is NE-layer code. Where NE's own files disagree (`TYPEDSIGNATURES` in `AGENTS.md`, `SIGNATURES` in `docstring-rules.md`), `AGENTS.md` wins: `TYPEDSIGNATURES`, and literal zeros (`max(0, a)`). Oceananigans' rules apply to anything sent upstream there. | — |
| D6 | Time stepper | Split Runge-Kutta (NE's default) on every level. It carries no tendencies from the previous step, so overwriting the parent's state between steps leaves nothing stale. | Adams-Bashforth, which would also need its stored tendencies reconciled after feedback. |

## 2. What already exists

**NumericalEarth 0.9.2, `NestedModels`** (one-way, written for Breeze atmospheres):
- `NestedModel(parent, child, exchanger, child_callbacks)`: an `AbstractModel` that forwards the
  model protocol to the child. Its `time_step!` steps the child, then advances the parent by
  the same `Δt`. There is no time refinement and no feedback.
- `Interpolated(source)`: a boundary-condition value that interpolates a `Field` (space) or a
  `FieldTimeSeries` (space and time) at the boundary face. `VectorComponent` rotates vectors.
- `parent_boundary_conditions`, `parent_forcings` (interior `Relaxation` toward a parent
  `FieldTimeSeries`), `blend_parent_terrain!`.
- `exchange_state!(exchanger, time)`: the hook a `NestedModel` calls to refresh parent-derived
  data. The Breeze extension's `StateExchanger` keeps a two-level time window and is the
  precedent for the ocean exchanger here.

**Oceananigans 0.117.2, open boundaries for `HydrostaticFreeSurfaceModel`:**
- `GravityWaveRadiationBoundaryCondition((U, η))` (Flather) on the barotropic transport, with
  optional `target_transport`; `SurfaceWaveRadiation` (Chapman) on `η`; both applied every
  barotropic substep.
- `TracerReservoir`: tracer open boundary with exterior value on inflow, interior on outflow.
- Normal and oblique radiation, `PerturbationAdvection`, `TidalHarmonics` with
  `tidal_boundary_conditions` and `tidal_forcing`.
- `TracerTendencyCallsite` callbacks, run after tracer tendencies in every Runge-Kutta stage.
- `RegriddedOperation` plus the ConservativeRegridding extension; `regrid!` in
  `Fields/regridding_fields.jl` handles one dimension at a time on matching grids.

**`~/NWEuropenShelf`** (the earlier prototype): a working one-way, barotropic, single-layer
chain with the recursive `advance!` ordering this plan adopts (parent steps, child subcycles
against time-interpolated parent data). It has no tracers, hand-rolled boundary arrays,
its own bathymetry regridding and interpolation stencils, and grids whose faces do not align
between levels.

**What is missing, and is the work of this plan:**
1. A nested model whose step order is parent first, child subcycled, then feedback, and which
   composes recursively.
2. An exchanger for ocean models: a two-time-level copy of the parent state around the child.
3. Aligned nest grids with consistent bathymetry.
4. Restriction of the child's state onto the parent.
5. Flux correction at the interface.
6. The three-level test case.

## 3. Design

### 3.1 Composition and step order

Nests are right-nested: `nest(shelf, nest(approaches, bay))`. A nest's child may itself be a
nest; its parent is always a single model. The opposite association cannot work, because the
bay has to subcycle inside each step of the approaches.

```
time_step!(nest, Δt)
    time_step!(parent, Δt)                    # parent to t + Δt
    exchange_state!(exchanger, t + Δt)        # parent window now brackets [t, t + Δt]
    for n in 1:refinement
        time_step!(child, Δt / refinement)    # recursion happens here
    end
    feed_back!(parent, child)                 # restriction, then update_state!(parent)
```

`NestedModel.time_step!` in NE is not dispatched on anything, so this needs either a new
method keyed on the exchanger type or a new model type. Default: a new type in this package,
`TwoWayNestedModel`, holding `parent`, `child`, `exchanger`, `refinement` and the parent's
clock. It is an `AbstractModel`, so `Simulation(nest; Δt)` works and `Δt` is the outermost
step. If the type ends up differing from `NestedModel` only in `time_step!`, milestone 9
proposes the upstream change instead of keeping both.

Two accessors keep levels generic (dispatch, no branching):
- `ocean_model(level)`: the `HydrostaticFreeSurfaceModel` of a level, whether the level is a
  bare model, an `EarthSystemModel`, or a nest (its parent).
- `outer_level(child)`: the level a parent exchanges with; the child itself, or a nest's parent.

### 3.2 Construction order

A child's boundary conditions read the exchanger, and the exchanger reads the parent, so
levels are built coarse to fine and assembled afterwards. The intended user script:

```julia
shelf_grid      = LatitudeLongitudeGrid(...)                                   # or any NE grid
approaches_grid = nested_grid(shelf_grid;      longitude, latitude, refinement = 5)
bay_grid        = nested_grid(approaches_grid; longitude, latitude, refinement = 5)
# bathymetry on each grid, reconciled fine to coarse (3.4), then immersed grids

shelf      = EarthSystemModel(...)                                             # any NE model
approaches = nested_ocean(shelf,      approaches_grid; kw...)                  # exchanger, boundary conditions, sponge, model
bay        = nested_ocean(approaches, bay_grid; kw...)

model      = TwoWayNestedModel(shelf, TwoWayNestedModel(approaches, bay))
simulation = Simulation(model; Δt, stop_time)
```

`nested_ocean` and how the exchanger travels with the child are a sketch; the names and
shapes are settled in milestone 2.

**A parent handed over as a finished model** has a fixed grid, so its bathymetry cannot be
reconciled with the child's. Nesting still works; the parent and child volumes differ inside
the nest, restriction is not exactly conservative, and the drift is measured and reported
rather than removed.

### 3.3 Downscale: boundary state of the child

The exchanger holds, for each exchanged parent field, a two-level in-memory `FieldTimeSeries`
on a window of the parent grid that covers the child plus a margin. `exchange_state!` shifts
level 2 to level 1 and copies the new parent state into level 2. Land points in the window
are filled from neighboring sea points (NE's `inpaint_mask!`) so interpolation near the coast
does not pull in zeros; this replaces the prototype's weight-renormalizing stencil.

Child boundary conditions are built with the existing `parent_boundary_conditions` and
`Interpolated`, which already interpolate a `FieldTimeSeries` in space and time:

| Child field | Condition | Exterior value |
|-------------|-----------|----------------|
| `U`, `V` | `GravityWaveRadiationBoundaryCondition` | parent `(U, η)` |
| `η` | `SurfaceWaveRadiation` (the default at open sides) | — |
| `u`, `v` normal | `NormalFlowBoundaryCondition` with a radiation scheme | parent `u`, `v` |
| `u`, `v` tangential | `ValueBoundaryCondition` | parent `u`, `v` |
| tracers | `TracerReservoir` | parent tracer |

A sponge of a few parent cells inside the child relaxes tracers and velocities toward the
parent through `parent_forcings`.

To verify in milestone 2, because reading the source did not settle them:
- that `Interpolated` works as an element of the Flather `(U, η)` tuple;
- that 3D normal velocities and the barotropic transport stay consistent at an open side with
  several layers (the prototype had one layer and a zero normal-velocity condition);
- how the Breeze `StateExchanger` rewrites its window's times when `Δt` changes and fills the
  first level at iteration 0, so the ocean exchanger follows the same pattern;
- whether `inpaint_mask!` (iterative nearest-neighbor fill) is cheap enough to run on 3D
  window fields every parent step, or whether the fill pattern should be computed once.

### 3.4 Grids and bathymetry

`nested_grid(parent_grid; longitude, latitude, refinement)` snaps the requested box to
parent faces and returns a `LatitudeLongitudeGrid` with `refinement` child cells per parent
cell in each horizontal direction and the parent's vertical faces, truncated at the child's
maximum depth. It also returns the parent index ranges the child covers, which the
restriction and the exchanger window both use.

Bathymetry:
- Each level regrids its own source data, through NE's `regrid_bathymetry`.
- Within the sponge the child's bottom is blended to the parent's (`blend_parent_terrain!`
  applied to bottom height), so the boundary state is consistent with the depth it came from.
- Inside the feedback region the parent's bottom height is replaced by the area mean of the
  child's, working from the finest level outward. Column volumes then match. With partial
  cells the volumes of individual levels near the bottom still differ, which is one source of
  the drift measured in milestones 4–5.

The 100 m level needs finer bathymetry than GEBCO. The EMODnet data the prototype used gets a
dataset type with the `DataWrangling` hooks in place of the prototype's `regrid_elevation`.

### 3.5 Upscale: state feedback

Applied once per parent step, after the child has reached the parent's time.

| Field | Restriction onto the parent | What it preserves |
|-------|-----------------------------|-------------------|
| tracers | `cᴾ = Σ cᶜ Vᶜ / Σ Vᶜ` over the wet child cells inside the parent cell | tracer content, constants |
| `u`, `v` | `uᴾ = Σ uᶜ Aᶜ / Σ Aᶜ` over the child faces lying on the parent face | volume transport through each parent face, hence the divergence |
| `η` | area-weighted mean of the child's `η` | volume |
| `U`, `V` | not restricted separately: recomputed from the restricted `u`, `v` by the model's own barotropic-mode computation | consistency of the two modes |
| `w` and other diagnostics | recomputed by `update_state!(parent)` | — |

Notes:
- Volumes and face areas come from the grid operators, so thickness-weighting under a z-star
  coordinate is automatic.
- Restricting `u`, `v` as face transports means only child faces that coincide with a parent
  face contribute, which is why aligned grids (D1) matter more for momentum than for tracers.
- The split-explicit free surface also carries a filtered (time-averaged) barotropic state
  that advects tracers. Whether the parent's filtered state has to be reset from the child's
  as well, so the parent's next tracer step sees transports consistent with the fed-back `η`,
  is settled in milestone 5 by the volume budget.
- The feedback region stops a few parent cells short of the child boundary by default (the
  sponge width), where the child is relaxed toward the parent anyway.
- For an `EarthSystemModel` parent, `update_state!` recomputes surface fluxes; if that cost
  matters, a narrower refresh is a later optimization.
- At initialization the child is restricted onto the parent once, so the two start consistent.
- Searched for an existing operator before proposing new kernels: `Fields/regridding_fields.jl`
  (`regrid!` needs matching grids in the other two dimensions), `RegriddedOperation` and the
  ConservativeRegridding extension (general, horizontal-area weights, CPU-built; not yet
  checked for immersed cells, volume weighting or face fields). Milestone 4 starts by testing
  whether the extension covers the tracer case; if it does, it is used and no kernel is added.

### 3.6 Flux correction (milestone 7)

Restriction alone does not conserve the composite content (parent outside the nest plus
child), because the parent and the child compute different fluxes through the interface. The
correction: accumulate the child's volume and tracer fluxes through each interface face over
its `refinement` steps, compare with the parent's fluxes through the same face over its one
step, and add the difference to the parent cell just outside the interface.

The tracer correction depends on the volume transports agreeing. For a uniform tracer the
correction reduces to the transport mismatch through the face, so correcting tracers without
matching transports would break uniformity. This is why momentum and sea-level feedback
(milestone 5) come before flux correction, and why the correction is applied to volume
(through `η`) and tracers together.

Oceananigans does not store tracer fluxes; they are computed inside the tendency kernels. They
can be recomputed on the interface faces from the advection scheme, the transport velocities
and the tracer, in a `TracerTendencyCallsite` callback. With the split Runge-Kutta stepper
every stage restarts from the state at the start of the step, so only the last stage's flux
should enter the final update, making the accumulation one flux evaluation per step rather
than a weighted sum over stages. That reading is from the source, not from a test; milestone
7 opens by checking it against a closed-box budget. Diffusive fluxes through the interface
are included the same way.

Flux correction requires the feedback region to reach the interface, so the rim of 3.5 is
revisited there.

## 4. Milestones

Each ends with checks that pass or fail. Idealized cases come before Pembrokeshire.

**0. Housekeeping.**
- `Project.toml` currently has `name = "NWEuropenShelf"` and that project's UUID; give the
  package its own name and UUID.
- `git init`; the restraint self-review reads `git diff`.
- Add `test/` with ParallelTestRunner, and a `CLAUDE.md` pointing at the governing rules.
- Check: `using OnlineNesting` loads and an empty test suite runs.

**1. Aligned nest grids.** `nested_grid` as in 3.4.
- Check: every `refinement`-th child face equals a parent face to round-off; the parent cell
  area equals the sum of its child cell areas.

**2. One-way nest of 3D models.** Exchanger, child boundary conditions, sponge, and the
`TwoWayNestedModel` step order with `feed_back!` doing nothing. Two idealized levels on a
flat-bottomed box with open sides.
- Check A, rest: a stratified ocean at rest stays at rest on both levels.
- Check B, wave: a barotropic gravity wave crosses the child without visible reflection, and
  the child's sea level matches the parent's in the overlap to a stated tolerance.
- Resolves the four items flagged in 3.3.

**3. Recursion.** Three idealized levels.
- Check: the clocks of all levels agree after every outer step, and check B holds at the
  innermost level.

**4. Tracer feedback.** Tracer row of 3.5.
- Check A: a uniform tracer stays uniform on every level.
- Check B: a passive blob advected out of the child and back in; report the composite tracer
  content over time.
- Check C: the parent's blob inside the nest matches the restricted child.

**5. Momentum, sea-level and barotropic feedback.** Remaining rows of 3.5.
- Check A: rest stays at rest with feedback on.
- Check B: a geostrophically balanced eddy resolved by the child but not the parent appears
  on the parent after feedback and stays balanced there (no burst of gravity waves).
- Check C: composite volume over time, reported; the wave test of milestone 2 holds or improves.
- Settles the filtered-state question in 3.5.

**6. Real bathymetry.** Dataset hook for EMODnet, blending in the sponge, parent bottom from
the child in the feedback region.
- Check: rest state stays at rest on the three real grids; column volumes of parent and child
  agree in the feedback region.

**7. Flux correction.** As in 3.6.
- Check A: the blob test of milestone 4 conserves composite tracer content to round-off, and
  composite volume likewise.
- Check B: a uniform tracer still stays uniform.

**8. The Pembrokeshire case.** Section 5.
- Check: tidal sea level at St Brides Bay and currents in Ramsey Sound agree between levels
  as in the prototype's `report_agreement`; a tracer released in the bay appears on the
  approaches and shelf grids; composite volume and tracer content are conserved; the run is
  compared against the same chain with feedback off.

**9. Upstreaming.** Propose to NumericalEarth whatever turned out general: the step order,
the ocean exchanger, `nested_grid`. Propose to Oceananigans anything that belongs there
(a native `Interpolated`, exposed interface fluxes).

## 5. Test case: shelf → approaches → bay

Boxes start from the prototype's and are snapped to parent faces by `nested_grid`.

| Level | Longitude | Latitude | Spacing | Columns (approx.) | Step |
|-------|-----------|----------|---------|-------------------|------|
| shelf | −10.5 to −2.5 | 49.0 to 54.5 | 2.5 km | 220 × 245 | Δt |
| approaches | −6.4 to −4.6 | 51.3 to 52.4 | 500 m | 250 × 245 | Δt / 5 |
| bay | −5.50 to −5.06 | 51.68 to 51.94 | 100 m | 305 × 290 | Δt / 25 |

- Refinement 5 in space and time at both nests.
- Vertical: one set of faces shared by all levels, surface-intensified, down to the shelf
  level's maximum depth; the inner levels use the upper part.
- Tracers: `T`, `S`, and one passive tracer released in the bay to show the upscale exchange.
- Outer boundary: Flather with tidal harmonics through Oceananigans' `TidalHarmonics` and
  `tidal_boundary_conditions`; `T`, `S` from GLORYS through NE's `DataWrangling`.
- Surface: prescribed atmosphere (ERA5 or JRA55), one `EarthSystemModel` per level (D3).
- Initial state: GLORYS on the shelf; inner levels interpolated from their parent.

Open points for this case:
- A tidal atlas for the outer boundary (TPXO or FES) needs a data source and license. Until
  then the prototype's uniform M2 oscillation stands in, which exercises the nesting but is
  not the real tide.
- "NW European shelf" in the prototype is the Celtic and Irish Seas box above. If the
  outermost level should be the whole shelf (roughly 20°W–13°E, 40–65°N), only the first row
  of the table changes.
- NE 0.9.2 has no ready-made constructor for an ocean with open boundaries from a dataset;
  the outer level's boundary conditions are assembled in the example script from the pieces
  above until that exists upstream.
- No wetting and drying: a minimum depth is imposed, as in the prototype.

## 6. Execution: sessions and agents

**State lives in the repository, not in a conversation.** Every session starts by reading
`CLAUDE.md` (which points at the governing rules), this plan, the status log below and
`git log`. Every session ends by committing the milestone and appending to the status log:
what was built, the check results with their numbers, and any departure from the plan.

**One fresh session per block of milestones.** A milestone's code and checks carry everything
the next one needs; the reading and dead ends that produced them do not.

| Session | Milestones | Depends on | Pause for review afterwards |
|---------|------------|------------|-----------------------------|
| A | 0, 1 | — | no |
| B | 2 | A | yes: the user-facing shape (`nested_ocean`, the exchanger) is fixed here |
| C | 3, 4 | B | no |
| D | 5 | C | yes: the hardest numerics |
| E | 6 | A only | no; can run alongside B–D in a git worktree |
| F | 7 | D, E | yes: conservation numbers |
| G | 8 | F, and the shelf extent and tide decisions | yes |
| H | 9 | G | each upstream proposal is confirmed before it is opened |

**Inside a session, the main session writes the numerics.** The step order, exchanger,
restriction kernels and flux correction are tightly coupled, and a brief to another agent
would lose the detail that makes them correct. Agents are used for three things:
- *Source reading.* A read-only search agent answers a specific question about Oceananigans
  or NumericalEarth internals with file and line references: the four items flagged in 3.3
  (session B), the filtered barotropic state (session D), the tracer flux functions and the
  Runge-Kutta stage weights (session F).
- *Independent review.* At the end of each milestone a fresh agent gets the diff and the
  restraint rules, not the author's conclusions, and reports what it would cut or question.
  In sessions D and F a second reviewer checks the numerics against section 3.
- *Independent side work.* Session E, and data sourcing for session G, run as a separate
  agent in a worktree while the main line continues.

Long simulations run in the background. No multi-agent workflow is planned; the work is
sequential enough that orchestration would cost more than it saves.

### Status log

- 2026-11-10: plan written and decisions settled. No code yet.
- 2026-11-10, session A (milestones 0 and 1):
  - Package renamed to `OnlineNesting` with its own UUID; `DocStringExtensions` added as a
    dependency for `TYPEDSIGNATURES`. Git repository initialized on `main`; `Manifest.toml` is
    ignored. `CLAUDE.md` and copies of NumericalEarth's `.claude/rules` added. Tests run
    through ParallelTestRunner with `Pkg.test()`.
  - `nested_grid(parent_grid; longitude, latitude, refinement, depth)` and
    `parent_indices(parent_grid, nest_grid)` in `src/nested_grid.jl`. Departure from 3.4:
    `nested_grid` returns only the grid, and the parent index ranges come from
    `parent_indices`, which derives them from the two grids.
  - Check passed, 56 of 56 assertions in `test/test_nested_grid.jl`: for regular, stretched,
    z-star and immersed parents at refinements 2 and 5, every `refinement`-th nest face matches
    a parent face (`≈`), the vertical faces are identical, and the nest cell areas inside each
    parent cell sum to its area (`≈`).
  - Not covered: GPU, and a nest whose requested box reaches the parent's edge (the extents
    snap to the nearest parent face with no check that a margin remains for the exchanger).

## 7. Risks

- **Open-boundary consistency with several layers.** The barotropic and baroclinic normal
  flow at an open side have to agree. Oceananigans has machinery for this, but it has not
  been exercised here beyond one layer. Milestone 2 is where this shows up.
- **Barotropic feedback under split-explicit stepping.** The parent's sea level, transports
  and filtered state are overwritten between steps. Getting these mutually consistent is the
  hardest part of two-way nesting in split-explicit models, and milestone 5 is the one most
  likely to take longer than it looks.
- **Feedback rim versus flux correction.** A rim keeps the feedback away from the sponge;
  flux correction wants the feedback to reach the interface. Milestone 7 has to reconcile them.
- **Cost of the time-interpolated boundary lookup.** `Interpolated` searches the parent
  window at every boundary point, stage and barotropic substep. If it shows up in profiles,
  precompute the boundary values once per child step.
- **GPU.** Everything is planned as kernels with on-device arrays, but land filling and the
  ConservativeRegridding weights are CPU-side today.
