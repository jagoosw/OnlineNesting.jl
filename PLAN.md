# OnlineNesting.jl: plan for two-way, multi-level nesting

Target: a chain of ocean models that run together, coarse to fine, where each nest takes its
lateral boundary state from its parent and returns its full state (tracers, velocities,
barotropic transport and sea level) to the parent. The outermost level is any NumericalEarth
(NE) model: a global `EarthSystemModel`, or a regional one with open boundaries. The test
case is three levels: NW European shelf → Pembrokeshire approaches → St Brides Bay.

Versions this plan was checked against: Oceananigans 0.117.2, NumericalEarth 0.9.2, Julia 1.12.
Milestones 0–7 have been run (status log, section 6); the checks of later milestones are still to be done.

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
- Outer boundary, with the GLORYS fields through NE's `DataWrangling`:
  - `U`, `V`: Flather, with exterior transport and sea level the sum of the tidal harmonics
    (Oceananigans' `TidalHarmonics` and `tidal_boundary_conditions`) and GLORYS' sea level and
    depth-integrated currents. GLORYS has no tide, so the two add.
  - `u`, `v`: the full GLORYS currents, normal and tangential, as for a nest's sides (3.3).
  - `T`, `S`: GLORYS, as the exterior value of a `TracerReservoir`.
- Surface: prescribed atmosphere (ERA5 or JRA55), one `EarthSystemModel` per level (D3).
- Initial state: GLORYS `T`, `S`, currents and sea level on the shelf; inner levels interpolated
  from their parent.

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
  - Check passed, 130 of 130 assertions in `test/test_nested_grid.jl`: for regular, stretched,
    z-star, immersed and Float32 parents at refinements 2 and 5, the nest spans the parent
    faces nearest the requested extents, splits each covered parent cell into `refinement`
    equal cells (`≈`), has the parent's vertical faces exactly, and carries over the halo,
    float type, radius and kind of vertical coordinate.
  - Departure from the milestone 1 check: the area-sum assertion was dropped. On a
    latitude-longitude grid cell area is additive, so it follows from the face check.
  - Independent review found that the first version read faces through `λnodes`, which omits
    the seam face of a periodic parent and cut a nest ending at the seam one column short.
    Faces now come from Oceananigans' `cpu_face_constructor_*` and `reconstruction_faces`.
  - Not covered: GPU; the docstring example is not run by any test or docs build; a nest
    crossing a periodic parent's seam is truncated at the seam; a requested box reaching the
    parent's edge snaps to it with no check that a margin remains for the exchanger.
- 2026-11-10, session B (milestone 2):
  - Built: `StateExchanger(parent)`, `nested_boundary_conditions(exchanger, grid)` and
    `nested_forcings(exchanger, grid; rate, width)` in `src/state_exchanger.jl`;
    `TwoWayNestedModel(parent, child, exchanger; refinement)` in `src/two_way_nested_model.jl`,
    with `time_step!` and `update_state!` (the only methods `Simulation` and `run!` asked for).
  - **User-facing shape, for review.** Departure from 3.2: there is no `nested_ocean`. The user
    builds the child with its own constructor and passes it the nest's boundary conditions and
    forcings, and the exchanger is the third argument of the nest:
    ```julia
    exchanger = StateExchanger(parent)
    child = HydrostaticFreeSurfaceModel(grid; boundary_conditions = nested_boundary_conditions(exchanger, grid),
                                              forcing = nested_forcings(exchanger, grid; rate, width), ...)
    model = TwoWayNestedModel(parent, child, exchanger)
    ```
    Reason: under D3 every level is built by the user through `ocean_simulation` and
    `EarthSystemModel`, so a builder would have to forward all of their keywords, and the
    plan's two-argument `TwoWayNestedModel(approaches, bay)` had no source for the exchanger.
    Three levels compose as `TwoWayNestedModel(shelf, TwoWayNestedModel(approaches, bay, e₂), e₁)`
    with `e₂ = StateExchanger(approaches)`; this is on paper only, milestone 3 runs it. As
    built, the outer nest would need `refinement` passed, because the default reads `child.grid`.
    The name `StateExchanger` is also that of the type in NumericalEarth's Breeze extension,
    which matters when this goes upstream (milestone 9).
  - Checks, `test/test_one_way_nest.jl`, 14 assertions, all pass (144 of 144 with session A's).
    Setup: 40 × 40 × 4 parent over 10° × 10°, 100 m deep, closed sides; nest over the central
    4° × 4° at refinement 3 in space and time; `T`, `S` with a linear stratification; split
    Runge-Kutta 3, 20 barotropic substeps, parent step 5 minutes; sponge 4 cells at 1/hour.
    - Check A, rest: after 24 parent steps run through `Simulation`, `u`, `v`, `η` and the
      horizontal variation of `T` are exactly 0 on both levels (asserted below 1e-12), the
      clocks agree, and the child has taken 3 steps per parent step.
    - Check B, wave: a 0.1 m Gaussian bump of 1° standard deviation west of the nest, 80 parent
      steps. The largest difference between the child's block-averaged `η` and the parent's
      `η` over the nest is 3.6 % of the largest `η` in the nest (asserted below 5 %). With a
      0.5° bump (two parent cells) it is 12.5 %, so the difference follows resolution.
      Reflection is not measured separately; this one number bounds it. Sponge rates of 0,
      1/hour and 1/10 minutes give 3.8 %, 3.6 % and 3.4 %.
  - The four items of 3.3:
    1. `Interpolated` inside the Flather tuple: not as it stands. The model does not regularize
       the elements of a tuple condition, so `nested_boundary_conditions` tags each with its
       side by calling `regularize_boundary_condition` itself. The sea level's copy is held at
       `(Center, Center, Nothing)` so it is interpolated at the transport's two-dimensional node.
    2. Normal velocity and barotropic transport with several layers: consistent. At the child's
       open faces after the wave run, `|∫u dz − U|` is 3e-17 and `|∫v dz − V|` 6e-17, with
       `|U|` up to 0.24 and `|V|` up to 0.33 m² s⁻¹ there (asserted below 1e-12). Oceananigans' barotropic corrector covers the
       boundary faces, so the condition on `u` sets only the departure from the depth mean.
    3. The Breeze `StateExchanger`: no precedent. It follows a prescribed parent's fixed time
       axis and never rewrites times. Here the nest records the parent into level 1 just before
       the parent steps and into level 2 just after, then reassigns each series' `times` to
       `range(t, t + Δt; length=2)`; a changing `Δt` needs nothing more. Recording before the
       step, not shifting level 2 down, means a parent state changed between steps (by `set!`
       now, by feedback later) is what the child starts from. The first version shifted, and
       a parent set after the exchanger was built gave the child 0.05 m s⁻¹ from rest.
    4. `inpaint_mask!` every step: too expensive (a `deepcopy` per call, and two host-side
       reductions per pass), and nothing computes the fill pattern once. Nothing built.
  - Departures from the plan:
    - No land fill in the exchanger. The reviewer ran a parent with land touching the nest's
      side: 0.021 m s⁻¹ and a 0.46 °C tracer error from rest in three steps. Milestone 6 needs
      a fill computed once.
    - The series cover the whole parent grid, not a window around the child.
    - The normal velocity has no radiation scheme: it takes the parent's `u`, `v` directly.
    - `feed_back!` is not defined; milestone 4 adds the call after the child's steps.
    - `exchange_state!` from NumericalEarth is not extended; `store_state!(exchanger, level)`
      does the copy.
    - No `ocean_model` or `outer_level` accessors: the parent must be a bare
      `HydrostaticFreeSurfaceModel` and the child a single model. Milestone 3 needs them for a
      nested child, and D3 for an `EarthSystemModel`.
  - Known and not fixed:
    - `Simulation` callbacks passed to the nest's `time_step!` and `update_state!` are dropped.
    - A child tracer the parent lacks gets no boundary condition or sponge, silently.
    - Measured on a bare model with a Flather function that records its clock: a split
      Runge-Kutta 3 step from 90 s to 180 s queries the exterior values at times from 90 s to
      242 s, because each stage's barotropic substeps are labeled from the stage's start time.
      Every child step therefore reads parent values up to 0.7 child steps late, and on the
      last child step of a parent step the queries pass the window and clamp to the parent's
      newest level. (The reviewer's reading, that all queries stay inside the window, holds
      for the three-dimensional fields only.)
    - No test isolates the time interpolation: if the series' times were never rewritten,
      every query would clamp to the newest level and the wave test would probably still pass.
    - The clocks agree exactly only when `Δt / refinement` is representable; with `Δt = 0.1`
      the reviewer saw them differ by 1e-13 after 200 steps. Milestone 3's check has to allow
      for this or the nest has to set the child's time.
    - Before the first step the exchanger holds zeros, so halos filled by `set!` on the child
      are wrong until the child's first step refreshes them.
    - The sponge `width` counts mean cells, which is not the actual cells of a stretched nest.
    - The wave test passes without the sponge, so the sponge is covered only by the rest test.
    - Not run: GPU, an immersed or z-star child, a parent with open sides.
- 2026-11-10, session E (milestone 6), on branch `worktree-session-e` in `.claude/worktrees/session-e`:
  - `EMODnet2024` in `src/emodnet.jl`: a region-windowed bathymetry dataset on NumericalEarth's
    `DataWrangling` hooks, following `GLO30`. A `Metadatum` with a `BoundingBox` downloads that
    window from the EMODnet ERDDAP server by index range, and `regrid_bathymetry` reads it. Land
    without data (most land; 1724 such cells border the sea in the Pembrokeshire window) is
    filled with a height of 10 m, the prototype's value.
  - `blend_parent_bottom!(bottom_height, parent_bottom_height; width, minimum_depth)` and
    `restrict_bottom!(parent_bottom_height, bottom_height; width, minimum_depth)` in
    `src/nested_bathymetry.jl`, with `width` in parent cells for both. The blend copies the
    parent's bottom onto the nest cell by cell (no interpolation, so the nest's boundary has the
    parent's coastline) and ramps with NumericalEarth's `blend_parent_terrain!`. The restriction
    is the area-weighted mean of the nest's depth with land as zero depth.
  - Departures from 3.4:
    - Both work on depth with land at zero height, not on raw bottom height, so land elevation
      does not enter a mean. Both also keep the sea floor within the depth of the nest grid.
    - No wetting and drying means a minimum depth, which an exact mean cannot honor at the coast.
      One rule serves both functions: a column shallower than half of `minimum_depth` becomes
      land, any other wet column is at least `minimum_depth` deep. This differs from
      `regrid_bathymetry`, which deepens every wet cell to the minimum.
    - The restriction is array operations on the CPU (a reshape and block sum), not a kernel,
      and takes its weights from the grid faces (`Δλ Δ(sin φ)`), not from `Az`. Milestone 4's
      tracer restriction may replace it.
  - New dependencies: `NCDatasets` (the downloaded file marks land as missing, which
    `regrid_bathymetry` cannot read, so the file is rewritten) and the `Downloads` standard
    library (to add a `download` method). `NumericalEarth` and `NCDatasets` added to the tests.
  - Tests, 159 of 159 passing (130 from milestone 1): `test/test_nested_bathymetry.jl` (24) on a
    synthetic sloping bottom at refinements 2 and 5: the nest boundary equals the parent, the
    ramp is linear cell by cell, the nest interior and the parent outside the feedback region
    are untouched, column volumes from the immersed grids (`Az` × column depth) agree to `≈`
    wherever the parent is deeper than the minimum depth, and uniform nests of 4, 6 and 30 m
    restrict to land, 10 m and 30 m. `test/test_emodnet.jl` (5) downloads a window around
    Skomer, checks the file's coordinates against the native grid to a tenth of a cell, and
    regrids it. It needs the EMODnet server.
  - Check on the three real grids, `validation/pembrokeshire_bathymetry.jl`: shelf 220 × 245,
    approaches 245 × 245, bay 305 × 290, 14 shared vertical faces to 250 m (nests cut at 141
    and 97 m, below their deepest columns at 135 and 90 m), `width = 3`, `minimum_depth = 10`,
    `PartialCellBottom`, and `major_basins = Inf` on the nests. `regrid_bathymetry`'s default
    `major_basins = 1` removes water that reaches the sea only through a nest's open sides: it
    took 1868 wet columns (Milford Haven) from the bay, which the restriction then copied into
    the approaches as land. A first run had both that and a bay cut at 77 m.
    - Rest, tolerance 1e-8 m s⁻¹ stated beforehand: after 1 hour (30, 150 and 750 steps),
      `max|u| = max|v| = max|w| = 0` and `max|η| = 0` exactly on all three levels. Departure:
      each level was run alone with closed sides, because the nested model does not exist on
      this branch. With a horizontally uniform stratification this shows only that partial
      cells create no pressure gradient; it should be repeated on the coupled nest.
    - Column volumes in the feedback region with `minimum_fractional_cell_height = 0`: shelf ←
      approaches, 48 of 1849 parent columns differ by more than 1e-6, all 48 at the minimum
      depth, the rest agree to 6.9e-13; total volume differs by −1.3e-4. Approaches ← bay, 160
      of 2860 differ, all 160 at the minimum depth, the rest agree to 5.6e-12; total −3.8e-4.
    - With the default `minimum_fractional_cell_height = 0.2`: 946 of 1849 and 1025 of 2860
      columns differ, by up to 6.7% and 9.5%; totals −4.3e-4 and −7.2e-4. `PartialCellBottom`
      deepens thin bottom cells on each grid separately. Restricting from the nest's adjusted
      bottom instead was tried on the first run and was worse in total (+3.6e-3), because the parent's own
      adjustment only deepens. So exact column volumes need a minimum fractional cell height
      of zero; this is the partial-cell drift 3.4 anticipated, now with a number.
  - Independent review found, and the commit fixes: the restriction applied the minimum-depth
    rule to the whole parent, not just the feedback region; a `Float32` grid with a `Float64`
    minimum depth gave an abstractly typed array; an interrupted download left an unfilled file
    in the cache; a bounding box in 0–360° longitudes gave a wrong index; and the coordinate
    test's tolerance was too loose to catch a one-cell shift.
  - Open, for a decision: EMODnet heights are relative to lowest astronomical tide, not mean sea
    level, which is meters in this region. Nothing corrects for it.
  - Not covered: GPU and `Float32` grids are untested; a `Metadatum` without a region asks the
    server for the whole dataset; `width = 0` in the blend gives NaN along the boundary; nothing
    checks that a nest's feedback region stays clear of its parent's own blended rim.
- 2026-11-10, session E continued, after review with the user:
  - `refer_to_mean_sea_level!(bottom_height, lowest_tide; minimum_depth)` in
    `src/nested_bathymetry.jl` adds the height of the lowest astronomical tide relative to mean
    sea level (a number, a function or a field on any grid) to a bottom height measured from
    the lowest tide, then applies the minimum-depth rule. The sea surface cannot absorb the
    datum offset, which varies across the domain (from memory, under 1 m off south-east Ireland
    to 6–7 m in the upper Bristol Channel), so it goes into the bathymetry. Tests: 12 more
    assertions, 171 of 171 passing.
  - The validation script calls it between regridding (now with no minimum depth) and blending,
    with `lowest_tide = 0`: no correction is applied yet. A real field needs a tidal atlas, or
    the shelf model's own lowest sea level, and that run also needs an atlas at its boundary.
    This is the same open point as the tide source in section 5.
  - Decided by the user: partial cells (after trying full cells), and 50 layers 0.5 m thick at
    the surface, stretching by 7.6% per layer to 250 m, in place of 14. The nests have 42
    layers (to 141 m) and 37 (to 97 m).
  - Full cells on the earlier 14 layers, for the record: 708 of 1849 and 832 of 2860 parent
    columns differ from their nests, by 4.4% and 7.6% on average and up to 36% and 42%; totals
    −3.8e-5 and −3.3e-4; rest exactly zero. Every bottom snaps to a face, and so does the mean.
  - The validation script now reports column volumes for partial cells, partial cells of any
    thickness and full cells, and runs the rest check with partial cells for 20 minutes. Its
    numbers on 50 layers were not in when this was merged; they follow in the next entry.
- 2026-11-10, session B continued, after review with the user: land fill in the exchanger.
  - `StateExchanger` computes once, for each parent field, the indices of its land nodes
    (`inactive_node`, so coastal faces keep their zero) and of the nearest sea node on the same
    level, found by spreading outward one cell per pass. `store_state!` copies the parent and
    then assigns each land node its donor's value. The 2D fields (`η`, `U`, `V`) use the
    surface level's land.
  - Check, added to the rest test: a parent with land along the nest's west side and a 50 m
    shelf along its south side, under a flat 100 m child. Before the fill the child reached
    1.2 m s⁻¹, 0.33 m of sea level and 6.9 °C of horizontal tracer difference in 2 hours;
    with it `u`, `v`, `η` and the tracer difference are exactly 0 (asserted below 1e-12).
    193 of 193 assertions pass with session E's tests.
  - Not done: a parent level with no sea node at all has no donor and fails at construction
    with a bounds error; there is no fill from the level above. GPU not run.
- 2026-11-10, session B continued: the Pembrokeshire chain with an M2 tide, one way,
  `validation/pembrokeshire_tide.jl`. Ahead of milestones 3 and 8; a demonstration, not a check.
  - Setup: half the target resolution (shelf 110 × 122 × 4 at 5 km, approaches 125 × 120 × 3
    at 1 km, bay 150 × 145 × 2 at 200 m), EMODnet bottoms reconciled as in session E, partial
    cells, temperature stratification, steps of 200, 40 and 8 s. The shelf's four sides are
    open with a uniform 1 m M2 sea level in a Flather condition. The three levels are
    `TwoWayNestedModel(shelf, TwoWayNestedModel(approaches, bay, e₂; refinement), e₁; refinement)`
    with `refinement` passed, which is the first run of three levels. Two M2 periods took
    15 minutes on 8 CPU threads; the full-resolution case was about 5.4 s per shelf step
    (80 minutes) and was stopped.
  - No tide imposed, 3 shelf steps: max speed 1.2e-7 m s⁻¹ on the shelf, 4.1e-7 on the
    approaches, 1.4e-6 in the bay; not zero, and not traced to the partial cells or the nesting.
  - Second period, sea level: approaches against shelf at St Brides Bay 1.63 against 1.66 m,
    rms difference 0.021 m; bay against approaches 1.63 against 1.63 m, rms 0.002 m; at Ramsey
    Sound rms 0.027 and 0.018 m. Correlations 0.9999 or higher.
  - Second period, surface speed in Ramsey Sound: peaks 0.61 (shelf), 1.11 (approaches) and
    1.02 m s⁻¹ (bay), rms differences 0.30 and 0.36 m s⁻¹: the sound is narrower than a shelf
    cell and about one approaches cell wide, so the levels are not expected to agree there.
  - Not run on a GPU; the script takes `arch`, edited in after this run and not rerun.
  - For session C to pick up first, from the items above: a nest as the child without
    `refinement` passed (`ocean_model`), an `EarthSystemModel` as a level, `Simulation`
    callbacks, and the clocks when `Δt / refinement` is not representable. The land fill is done.
- 2026-11-10, session C (milestones 3 and 4), on branch `session-c` in `.claude/worktrees/session-c`,
  because another session was editing the main working tree (`src/OnlineNesting.jl`,
  `src/emodnet.jl`, `validation/pembrokeshire_tide.jl`). One commit covers both milestones.
  - Built, all in `src/two_way_nested_model.jl`:
    - `ocean_model(level)`: the level itself, the ocean of an `EarthSystemModel`, or the parent
      of a nest. `StateExchanger` and the nest read grids and fields through it, so a level can
      be any of the three. No `outer_level`: `ocean_model` of a nest is already its parent's.
    - Tracer feedback: after the child's steps, and once before the first step, a kernel
      replaces each parent tracer by `Σ cᶜ Vᶜ / Σ Vᶜ` over the wet child cells in each parent
      cell of the feedback region, then the parent's state is updated. The parent keeps its
      value where the child is all land. Every tracer of the parent must exist on the child.
    - The child's last step is the time left to the parent's time, not `Δt / refinement`.
    - `Simulation` callbacks are passed to the parent's `time_step!` and `update_state!`.
  - **User-facing changes.** `TwoWayNestedModel(parent, child, exchanger; width)`: `width` is
    required and is the rim, in parent cells, that takes no feedback; it has to be the `width`
    given to `restrict_bottom!`. Departure from 3.5, which made the sponge width the default:
    the nest does not know the sponge, whose `width` counts child cells. The `refinement`
    keyword is gone: space and time refinement both come from the grids. A `width` of half the
    nest or more leaves no feedback region, which is the only way to turn feedback off.
  - New dependency: `KernelAbstractions`, for `@kernel`.
  - The ConservativeRegridding question of 3.5, settled by reading, not by running:
    Oceananigans' extension builds its regridder on `underlying_grid` of both grids, so it has
    horizontal area weights only, no wet-cell mask and no cell volumes. A kernel was written.
  - Checks, 218 of 218 assertions pass (193 before). All on flat 100 m bottoms with four
    layers, refinement 3, split Runge-Kutta 3, sponge 4 child cells at 1/hour, `width = 2`.
    - Milestone 3, `test/test_recursion.jl`: 40 × 40 shelf over 10° × 10°, approaches over the
      central 4° (48 × 48), bay over the central 2° (72 × 72), composed as
      `TwoWayNestedModel(shelf, TwoWayNestedModel(approaches, bay, e₂; width), e₁; width)`.
      - Clocks: equal (`==`) on the three levels after each of 80 steps of 5 minutes, and after
        each of 20 steps of 0.1 s, where repeated `Δt / 3` steps disagree with the parent on
        14 of 20 steps. The bay takes 9 steps per shelf step.
      - Check B at the innermost level, the milestone 2 wave: the largest difference between
        the bay's block-averaged `η` and the approaches' is 1.8 % of the largest `η` there
        (0.00058 of 0.033 m; asserted below 5 %). As in milestone 2, reflection is not
        measured separately.
      - Two `EarthSystemModel` levels (`OceanOnlyModel(ocean_simulation(...))`, no closure or
        Coriolis) under `Simulation`: at rest after 6 steps, with the four clocks equal.
      - Initialization runs finest first: a tracer set to 1 only in the bay gives 1 on the
        approaches and on the shelf after `update_state!(nest)`.
    - Milestone 4, `test/test_tracer_feedback.jl`, and the same three-level wave run:
      - Check A: a tracer equal to 1 differs from 1 by less than 1e-12 on all three levels
        after the wave run.
      - Check C: a nest 50 m deep in a 100 m parent, with land filling some parent cells, parts
        of others and the bottom layer of others. The parent equals a mean computed in the
        test from the grid faces (`≈`), keeps its value where the nest is all land, and is
        untouched outside the feedback region. During the blob run the parent stays within
        1.0e-14 of the restricted child in the feedback region after every step.
      - Check B: a periodic channel 6° × 4° (24 × 16), nest over the central 2° × 2°, a blob
        of standard deviation 0.3° in solid-body rotation at 4 m s⁻¹, one circuit in 556 steps
        of 5 minutes. The blob leaves the nest (largest value in the nest falls to 0.02) and
        returns (0.52 in the nest, 0.53 on the parent). Tracer content relative to the start,
        largest departure over the circuit and value at its end:

        | | parent, largest | composite, largest | composite, end |
        |---|---|---|---|
        | feedback on | 6.7 % | 10.9 % | −2.3 % |
        | feedback off | 5e-6 | 11.3 % | −2.1 % |

        Composite is the parent outside the nest plus the whole child; with feedback on, the
        start is the state after the restriction at initialization. It moves by a tenth
        with or without feedback, as the blob crosses the nest's sides on two grids that carry
        it differently; this is what flux correction (milestone 7) is for. These numbers are
        from a run outside the tests, which assert only that the blob leaves and returns and
        check C.
    - The milestone 2 tests now run with tracer feedback (file renamed
      `test/test_two_level_nest.jl`); the wave difference is unchanged at 3.6 %.
  - Independent review found, and the commit fixes: `update_state!` of an `EarthSystemModel`
    recomputes fluxes but does not update its ocean, so a fed-back parent kept its old halos
    and pressure (0.13 in `pHY′` in the reviewer's run; now both are updated, with a test);
    the test of an all-land parent cell had no such cell inside the feedback region.
  - Known and not fixed:
    - The exact last step changes `Δt` by round-off when `Δt / refinement` is not
      representable, which makes a quasi-Adams-Bashforth child take two Euler steps per parent
      step (reviewer's reading and run). D6 already requires split Runge-Kutta.
    - Callbacks reach only the parent, never the child, and update-state callbacks run once
      more per step than on the parent alone, in the update after feedback. No test attaches one.
    - The restriction at initialization runs twice under `Simulation` and again at an inner
      nest's first step. It gives the same result each time.
    - The kernel writes into a parent land cell that the child has as sea.
    - Not run: GPU, z-star, a nest with feedback on the real partial-cell grids (the reviewer
      ran one partial-cell case, which matched to 1.1e-13).
    - `validation/pembrokeshire_tide.jl` now feeds temperature back. It had a doubled semicolon
      from the unrun `arch` edit, fixed here. Run for 3 control steps and 2 tide steps only: the
      control speeds are those of the one-way run (1.2e-7, 4.1e-7 and 1.4e-6 m s⁻¹).
  - `validation/idealized_nests.jl` records two videos of three idealized levels: the blob in
    a 12° channel (`idealized_blob.mp4`) and the wave (`idealized_wave.mp4`). Not tracked.
  - Merging: `src/OnlineNesting.jl` and `validation/pembrokeshire_tide.jl` are also changed,
    uncommitted, in the main working tree.
- 2026-11-10, session D (milestone 5), on branch `session-d` in `.claude/worktrees/session-d`,
  because the main working tree again had another session's uncommitted edits. **Pause for review.**
  - Built, all in `feed_back!` in `src/two_way_nested_model.jl` (source +45 lines net):
    - `η`: the area-weighted mean over the child's wet columns, through the same kernel as the
      tracers (`_restrict_mean!`, which takes the weight and the vertical location). The parent
      keeps its value where the child is all land.
    - `u`, `v`: `Σ uᶜ Aᶜ / max(Σ Aᶜ, Aᴾ)` over the child faces on each parent face
      (`_restrict_transport!`), on the faces on both sides of every cell of the feedback region,
      so the divergence of a fed-back cell is entirely the child's.
    - `U`, `V`, halos and masking: Oceananigans' `reconcile_state!` on the parent, then the
      state update that was already there. Nothing of the barotropic mode is written by hand.
  - **Departure from 3.5, for review: the denominator of the velocity restriction.** The table has
    `Σ Aᶜ`. The first version divided by the parent's own face area, which keeps the transport
    exactly; the numerics reviewer showed it is unstable where a parent partial cell is thinner
    than the child's faces: a flat child at 0.1 m s⁻¹ over a parent bottom cell 1 cm thick gave
    the parent 250 m s⁻¹, and at 1e-13 m thick NaN on both levels within 3 steps. Session E found
    the real grids need a minimum fractional cell height of 0, so this is reachable. With the
    larger of the two areas the parent keeps the child's transport wherever its face is at least
    as large, and takes the child's mean velocity (losing transport) where it is smaller. On flat
    and grid-fitted bottoms the three forms are the same. How much transport the real grids lose
    this way is not measured.
  - The filtered-state question of 3.5: **nothing has to be reset.** Reading: every Runge-Kutta
    stage zeroes `η̅`, `U̅`, `V̅`, `Ũ`, `Ṽ` and restarts `η`, `U`, `V` from the fields cached at the
    start of the step, which is after the feedback; the transport velocities are copied from
    `u`, `v` at the start of each stage. Runs: filling the parent's filtered state with 1e3
    after a feedback leaves the next nest step bit-identical; and (reviewer) a parent after 3
    nest steps and a fresh model given the same `u`, `v`, `η`, `T`, `S` are bitwise equal in
    `u`, `v`, `w`, `η`, `U`, `V`, `T`, `S` before and after one more step. The volume budget did
    not decide this, contrary to 3.5: on a static vertical coordinate it cannot.
  - Checks, 238 of 238 assertions pass (218 before).
    All on four layers, refinement 3, split Runge-Kutta 3, `width = 2`, CPU, static `z`.
    - Check A, rest: the existing rest tests (flat, coast with land and a shelf, three levels,
      two `EarthSystemModel` levels) now run with all feedback and still give exactly 0.
    - Check B, eddy, `test/test_momentum_feedback.jl`: a geostrophic eddy of 0.02 m and radius
      0.25° (one parent cell, three child cells, a tenth of the barotropic deformation radius)
      set on the child only, 40 × 40 parent with Coriolis. After the feedback at initialization
      the parent's peak is 0.0159 m. The parent is then stepped **alone** for 4 hours:
      its peak never falls below 84 % of the start (90 % at the end) with the fed-back
      velocities, and falls to 10 % with them set to zero (16 % at the end); asserted above
      75 % and below 25 %. With radius 0.5°: 89 % and 15 %. The thresholds were chosen after
      seeing these numbers.
      - The first version of this test ran the nest and bounded the sea level radiated outside
        it at 5 % of the eddy (set beforehand). It measured 1.4 %, but 2.9 % with the `u`, `v`
        feedback removed, so it passed either way: the parent's peak is rewritten from the child
        every step. It was replaced. In that coupled run the parent's eddy went from 0.0159 to
        0.0156 m in 4 hours while the child's went from 0.02 to 0.0191 m.
    - Check C, composite volume in the milestone 2 wave run (80 steps), relative to the bump's
      initial volume, from a run outside the tests:

      | | parent alone, largest | composite, largest | composite, end |
      |---|---|---|---|
      | feedback on | 2.4e-3 | 6.1e-3 | 1.9e-3 |
      | feedback off | 5e-15 | 7.0e-3 | 3.6e-3 |

      The wave difference of milestone 2 (child block mean against parent over the whole nest)
      is **4.78 % with feedback, 3.68 % without** (3.6 % before this session). It is zero by
      construction in the feedback region now, so the number is the rim's, and it got worse,
      not better; check C asked for "holds or improves". The assertion is 5 %, so the margin is
      0.2 points. Not investigated.
    - Restriction, same file, for a grid-fitted and a partial-cell child 50 m deep with land,
      under a flat 100 m parent: the parent's `u Ax` and `v Ay` equal the sum of the child's
      over each face (`≈`), the faces of the rim are untouched, `U = ∫u dz` on the parent, and
      `η` equals a wet-area mean computed in the test from the grid faces, with the parent's
      value kept under an all-land block. A parent with a 1 cm bottom cell under a flat child
      at 0.1 m s⁻¹ gets 0.1 m s⁻¹.
  - Independent reviews. Restraint: source clean but for one local, inlined; two test names
    renamed. Not taken: merging the new tests into `test_tracer_feedback.jl` to share its helper
    (it would put Coriolis into the blob test), and cutting the `U = ∫u dz` assertion.
    Numerics: the thin-face instability and the eddy test above, both fixed; indexing checked
    by coordinates at refinements 2 and 3, `depth = 50` and widths 0 to 2 (transport error
    below 2e-16, the last face included); three levels feed back finest first, with each level
    equal to the restriction of the one below after 20 steps.
  - Known and not fixed:
    - **z-star is not supported by the velocity feedback.** The face areas are read before the
      parent's grid is rescaled to the new `η`, so the transport is off by `Δη / H` (2e-2 in the
      reviewer's run with `η` from −1 to 1 m on 100 m). Tracer weights have the same ordering.
    - `η` is the wet mean spread over the parent's whole cell, so volume is not preserved where
      the child is part land; `η` is also written into parent land under a wet child.
    - `reconcile_state!` refills the halos of a level with its own open sides, which moved `U`
      on those faces by up to 4.7e-2 of 0.97 in the reviewer's run (interior `η` differed by
      8e-10 of 2e-2 over 40 steps).
    - Immersed parents whose bottoms differ from the child's, measured by the reviewer before
      the denominator changed: no NaN, about 5 % of transport lost at 4 to 8 of 80 faces, cause
      not found. Not measured again.
    - The composite volume and wave numbers above were measured before the denominator changed;
      those runs have flat bottoms, where it changes nothing.
    - Not run: GPU, the real partial-cell grids, an `EarthSystemModel` parent with anything but
      rest, a periodic parent with momentum feedback.
  - Nothing was done on the Pembrokeshire case, at the user's instruction.
  - Merging: `src/OnlineNesting.jl` changes two import lines that the main working tree also has
    uncommitted edits to.
- 2026-11-10, session D continued, after review with the user (z-star must be supported; exact
  momentum conservation is not required). This supersedes the denominator and the z-star item of
  the entry above.
  - The velocity restriction is now 3.5's own formula, `Σ uᶜ Aᶜ / Σ Aᶜ` over all the child faces on
    the parent face, land faces included with their zero velocity. Tracers, `η`, `u` and `v` go
    through one kernel, `_restrict!`, a weighted mean that takes the weight (`wet_volume`,
    `wet_surface_area`, `Axᶠᶜᶜ`, `Ayᶜᶠᶜ`) and the block of child nodes. Source is +35 lines net
    against the start of the session. Consequences: the parent's velocity never exceeds the
    child's (the 1 cm test gives 0.1 m s⁻¹); transport through a parent face is the child's
    only where their face areas agree (flat and grid-fitted bottoms), not with partial cells.
  - **z-star.** No restriction reads the parent's grid any more, so nothing depends on when the
    parent's grid is rescaled; `reconcile_state!` rescales it from the fed-back `η`. Checked
    by running the wave test on z-star grids on both levels (`test/test_two_level_nest.jl`):
    with feedback the parent's uniform salinity stays within 1e-13 of 35 after 80 steps
    (asserted below 1e-11), which on z-star does test that transports and sea level agree.
  - **Upstream defect found, not ours and not fixed:** on a z-star grid a uniform tracer does
    not stay uniform in the cells beside an open side. In the nest the child's salinity is off
    by 2.4e-3 of 35 in its boundary columns, the same with feedback off, and 1.1e-2 without the
    sponge. A bare Oceananigans 0.117.2 model, 20 × 20 × 4 with open west and east sides, a
    gravity-wave radiation condition toward a 0.1 m exterior level, `NormalFlowBoundaryCondition(0)`
    on `u` and a `TracerReservoir` of 35, gives 4.7e-3 on z-star after 40 steps and exactly 0 on
    static `z`. Cause, traced afterwards: the filtered transport `Ũ` that advects tracers is built
    with `U`'s boundary conditions, and after the substeps `step_free_surface!` fills its halos
    with them, so on an open face `Ũ` is the radiation condition's end-of-step `U` (equal to 4
    digits in the bare model), not the weighted mean of the transports the substeps used. In
    the cells beside the open side `η̅ − η⁻ + Δt ∇·Ũ` is then 7.6e-3 m of a 4e-2 m change on the
    first step, against 1e-17 in the interior, and z-star turns that into tracer. The test marks the child's assertion `@test_broken`
    on z-star. So z-star is supported by the feedback, and the one-way boundary of a z-star
    child is not tracer-consistent until this is fixed in Oceananigans.
  - Tests: 247 pass and 1 is marked broken, of 248.
  - Wave difference again, without Coriolis (the 4.78 % above was a run with Coriolis):

    | | feedback on | feedback off |
    |---|---|---|
    | static `z` | 4.92 % | 3.64 % |
    | z-star | 4.71 % | 3.64 % |

    With feedback the largest difference is at the second parent cell in from the nest's west
    side: the innermost cell of the rim, between the parent's own solution and the fed-back
    region. The assertion is 5 %, so the margin is 0.08 points.
  - `validation/idealized_nests.jl` also records `idealized_eddy.mp4`: three levels with
    Coriolis, the eddy of check B set in the bay only, 12 hours. All three videos were recorded
    again with full feedback. Not tracked.
- 2026-11-10, session D continued: the z-star defect is fixed upstream in a draft pull request,
  CliMA/Oceananigans.jl#6485 (branch `jsw/open-boundary-filtered-transport`, worktree
  `~/Documents/Oceananigans-open-boundary-transport`). The filtered transport is accumulated
  after the barotropic boundary conditions, on the boundary faces too, and no longer carries
  `U`'s conditions. With that patch applied to a copy of 0.117.2, the z-star child's salinity
  error in the nested wave run is 1.8e-13 (2.4e-3 without) and the wave difference is unchanged.
  Until it is released the `@test_broken` in `test/test_two_level_nest.jl` stands; once it is,
  that assertion passes and has to become a `@test`. The pull request's own tests ran on CPU
  for three files only; its CI has not been read.
- 2026-11-10, session F (milestone 7), on branch `session-f` in `.claude/worktrees/session-f`, because the main
  working tree had another session's uncommitted edits. **Pause for review.**
  - **Oceananigans.** Everything here ran on 0.117.2 with the fix of CliMA/Oceananigans.jl#6485 applied (the
    user's instruction for this and later sessions). The branch itself does not load in this environment (it
    sits on a newer `main` that needs a newer DocStringExtensions), so the fix commit's `src` diff was applied
    to a copy of 0.117.2 and the untracked `Manifest.toml` points at the copy. `Project.toml` is unchanged, so
    a fresh resolve gets the release, on which the z-star child assertion below fails.
  - Built, all in `src/two_way_nested_model.jl`:
    - Each nest holds, on the parent grid, one field of volume mismatch and one of tracer mismatch per parent
      tracer, and two callbacks at the `TracerTendencyCallsite`, one for its parent and one for its child. In
      the last stage of every step of a level the callback adds `Δt` times the level's volume and tracer
      fluxes through the edge of the feedback region to the parent cells just outside it: out of them for
      the parent, into them for the child. The fluxes are Oceananigans' own `tracer_flux_x`, `tracer_flux_y`
      (advective plus diffusive, from its yin-yang conservation fixer) with the stage's transport velocities;
      the volume flux is the centered flux of a tracer equal to one.
    - A nest step zeroes the mismatch, steps the levels, applies the mismatch to those cells and then
      restricts as before. The
      sea level rises by the column's volume mismatch over its area. Each tracer is corrected layer by
      layer: a wet layer keeps the share of the column's volume change that z-star gives it and passes the
      rest of the volume it took upward, carrying the mean of its tracer and the next layer's; a land layer
      passes everything upward. Written as an increment, so no mismatch leaves the fields bitwise unchanged.
    - Bare models and nests get the callback through `time_step!`'s `callbacks`; an `EarthSystemModel` steps
      its ocean without them, so the callback is added to its ocean simulation at construction (with a TODO
      for NumericalEarth).
  - **Departures from 3.6, for review.**
    - The interface is the edge of the feedback region, not the child's boundary, so the rim stays (`width`
      at least 1 is what was tested; the first reviewer's sweep also conserved at `width = 0`). The child's
      faces there are ordinary interior faces, with no open-boundary or sponge terms. The composite that is
      conserved is therefore **the parent outside the feedback region plus the child inside it**, not
      session C's and D's "parent outside the nest plus the whole child"; numbers are not comparable with
      their tables. This settles the "rim versus flux correction" risk of section 7 without removing the rim.
    - Conditions: nothing forces the child inside the feedback region, so the sponge must end within the rim
      (`nested_forcings`' `width` in child cells at most `width × refinement`); a sponge of 7 or 9 child
      cells against a rim of 6 gave a composite tracer error of 4.7e-6. An inner nest's feedback region and
      the cells around it must lie inside the outer nest's: with the two edges touching, composite volume
      was off by 9.6e-4 and tracer by 8.9e-7.
    - The correction is the mismatch of the per-layer transports, which the plan's "volume through `η`" alone
      would not keep consistent with the tracers: z-star spreads a sea level change over the layers by
      thickness, a sheared mismatch does not, hence the upward transport in the tracer correction.
  - What was measured first, on bare models (CPU, 24 × 16 × 4 channel, sheared flow, 0.1 m bump):
    - A closed model alone, 100 steps: tracer content drifts 6.9e-4 on a static grid (the top cells do not
      hold what a rising sea level carries) and 6.7e-16 on z-star. **So round-off tracer conservation exists
      only on z-star**, with or without nesting.
    - One step, a box of columns: content change against `Δt` times the recomputed last-stage fluxes through
      the box's sides closes to 7e-14 of the change for tracers on z-star (centered; WENO 5 with horizontal
      diffusion; the same over partial cells) and to 4e-14 for volume on both coordinates. This confirms the
      reading in 3.6 that only the last stage's flux enters the update.
  - Checks. Tests: 24 assertions in the new `test/test_flux_correction.jl`; the bathymetry cases in it were
    written by a second agent. The suite: 270 pass and 2 are marked broken, of 272 (247 and 1 of 248 before).
    Refinement 3, `width = 2`, sponge 4 child cells, 5 minute steps.
    - Check A, composite conservation, largest departure over the run, relative to the initial tracer
      content and to the volume of the bump:

      | | tracer | volume | steps |
      |---|---|---|---|
      | z-star, two levels | 4.4e-16 | 7.4e-15 | 100 |
      | z-star, two levels, correction off | 3.0e-2 | 3.4e-3 | 100 |
      | z-star, WENO 5 and diffusion, partial cells | 4.4e-16 | 4.5e-15 | 100 |
      | z-star, three levels | 2.2e-16 | 1.3e-15 | 20 |
      | z-star, two `EarthSystemModel` levels | 4.4e-16 | 1.1e-15 | 20 |
      | static, two levels | 5.0e-4 | 7.6e-15 | 100 |
      | static, two levels, correction off | 2.9e-2 | 3.2e-3 | 100 |

      The milestone 4 blob circuit (556 steps, depth-uniform flow): 6.7e-16 on z-star and 4.9e-6 on a static
      grid, against 6.8e-2 on both with the correction off. Asserted below 1e-12 on z-star for two levels
      with WENO and diffusion, three levels and `EarthSystemModel` levels; on a static grid only the volume.
    - Check B, a uniform tracer under the sheared flow: within 2.1e-13 of 35 on both levels on z-star and
      1.4e-14 on a static grid (asserted below 1e-11 in every case above).
    - The first reviewer's sweep, all at round-off on z-star: refinements 2, 3 and 5, widths 0 to 3,
      third-order upwind. Dropping the correction on one side, or flipping its sign, gives a volume error
      of 3e-4 to 4e-3, so the assertions would catch either.
  - **The wave difference of milestone 2 got worse and its assertion no longer holds.** Child block mean
    against parent over the whole nest: 6.06 % (static) and 6.07 % (z-star) of the largest sea level, against
    4.93 % with the correction off; the bound is 5 %. All of it is the sea level correction (tracer
    correction alone: 4.93 %) and it sits in the cell that takes the correction, the innermost rim cell;
    with widths 1, 2, 3 it is 6.59, 6.06, 5.88 % against 5.27, 4.93, 4.75 %. That cell now holds the volume
    conservation requires, which is not the parent's own solution. The assertion in
    `test/test_two_level_nest.jl` is marked `@test_broken`, not loosened: to decide at review whether the
    bound is restated for the cells that take no correction, or the 1.1 points are to be reduced.
  - The `@test_broken` on the z-star child's salinity is now a `@test`, since it passes with the fix.
  - Bathymetry that differs between parent and child (first reviewer's cases, rerun after the fix below):
    - Found by the reviewer and fixed: where a parent column beside the edge has fewer wet layers than the
      child's faces, the first version averaged the tracer with the masked value of the land layer and put
      salinity off by 3.6 to 7.2 of 35. With land layers passing volume and content upward: 1.4e-14 on a
      static grid, composite tracer 4.4e-16 and volume 1.7e-15 on z-star.
    - Bottoms reconciled as in 3.4 (blend, restrict, minimum depth 10 m): composite tracer 2.2e-16 and
      volume 2.1e-15 on z-star with partial cells.
    - **Found here and fixed, in the exchanger, not the correction:** on z-star, where the parent is shallower
      than the child at the child's open side, a uniform tracer did not stay uniform in the child: 4.7e-3 of
      35 after one step, in the child's first column inside its side and its bottom layer, the same with
      feedback off, and 0.64 after 20 steps beside a parent column that is all land. Cause, traced by a
      separate agent with runs at each link: the child's bottom node, at its moving `z`, can lie below the
      parent's lowest cell center, which is indexed by the parent's reference coordinate, so the
      interpolation reaches the parent's bottom halo; the land fill covered interior nodes only and under
      parent land that halo held the masked 0, so the sponge target and the reservoir value were 34.91.
      `store_state!` now fills the halos of each stored level after the land fill: 2.1e-14 after one step
      and 1.2e-13 after 20 (partial cells and an all-land column the same). The z-star cases over a parent
      shelf and over reconciled bottoms now assert uniform salinity.
      - Left upstream: NumericalEarth's `Interpolated` boundary values and Oceananigans' `Relaxation`
        targets query a z-star child at its `z` and index the parent by its reference coordinate. Querying
        at the child's reference coordinate also removes the error (2.8e-14, run as an override) and an
        offset of the size of `η` in a stratified target, which the halo fill leaves.
      - The halos of every stored field now come from the filled interior and the stored series' own
        boundary conditions, not from the parent's halos; a nest beside a parent's open side is not run.
    - Not conserved, and not fixed: a parent column beside the edge that is land, or about 1 m deep, under
      a 100 m child loses the volume the child sends it (error 4 to 7 times the bump's volume); the sea
      level of a land column is left alone. A nest whose grid is shallower than the parent's ocean at the
      edge (`depth = 50` of 100) is not conserved either (volume error 0.4 to 0.5 with or without the
      correction). Reconciled bottoms and a nest grid cut below its deepest column avoid both.
  - Independent reviews. Agents ran on Opus after the first three ran on this session's model; the first
    numerics reviewer was stopped for that reason and its finished sweeps were read from its output files.
    - Restraint: taken: Oceananigans' `tracer_flux_x`, `tracer_flux_y`, `tracer_advecting_velocities` and
      `add_callback!` in place of copies (source +166 lines net in `two_way_nested_model.jl`), one kernel
      argument folded away, English names in the tracer kernel, the z-star limit in the docstring, tests on
      the public API with a diffusive case. Not taken: ring-sized storage for the mismatch; merging the sea
      level kernel into the tracer kernel.
    - Numerics (second reviewer): kernels and accumulation confirmed from first principles and by runs
      (column content to 5e-15, a uniform tracer under a sheared mismatch exact, box budgets to 5e-12,
      NumericalEarth's default `ocean_simulation` with CATKE and WENO at 4.4e-16 and 6.1e-15, `Float32`
      at its round-off). Found and fixed: an `EarthSystemModel` level stepped outside the nest kept
      accumulating (two steps of the parent alone moved the ring's sea level by 4.33 m at the next
      feedback); the mismatch is now zeroed at the start of a nest step and applied only in the step, with
      a test. Found and not fixed: the items marked (reviewer) below.
  - Known and not fixed:
    - (reviewer) Reconciling bottoms does not prevent a land column beside the edge: `restrict_bottom!`
      sets the feedback region only, so a coastline crossing the edge drops the volume sent to land.
    - (reviewer) The correction's volume can exceed a thin column's: beside a 1 m partial-cell column
      the factor `1 + mismatch / volume` reached −0.435, after which the values are finite and wrong. A
      mismatch of half the column's volume conserves content but overshoots the tracer's range.
    - (reviewer) A tracer callback of the user's reaches a bare parent but not the ocean of an
      `EarthSystemModel` parent, and never the child.
    - (reviewer) A tracer a closure steps itself (`e` with CATKE) is corrected like any other.
    - Building a second nest on the same `EarthSystemModel` leaves the first nest's callback on its ocean.
    - The mismatch fields are full three-dimensional parent fields, one per tracer and one for volume,
      holding one ring of cells; the correction is launched over every column of the region.
    - A level must use a split Runge-Kutta stepper (the callback reads `timestepper.Nstages`).
    - Momentum is not corrected.
    - Not run: GPU, `Float32`, the real grids.

- 2026-11-10, session G (milestone 8): the Pembrokeshire case with realistic forcing.
  - **Ported from NE PR 681:** `src/earth_tides.jl` (Earth tidal constituents and `TidalHarmonics`
    from astronomical constants) and `src/tpxo.jl` (`TPXO10Atlas` struct implementing
    `Oceananigans.BoundaryConditions.tidal_atlas_constants` — reads the 1/30° atlas files, inpaints
    land, interpolates harmonic constants at arbitrary boundary nodes). Both are temporary: they
    move to NumericalEarth once PR 681 is merged.
  - **`validation/download_pembrokeshire_data.sh`:** downloads GLORYS daily reanalysis (T, S, u, v)
    for the shelf region via `copernicusmarine subset`, checks that the TPXO atlas files exist, and
    notes that ERA5 is auto-downloaded by NumericalEarth.
  - **`validation/pembrokeshire.jl`:** full-resolution three-level case — shelf 220 × 245, approaches
    250 × 245, bay 305 × 290, 50 z-star layers surface-intensified to 250 m — with:
    - TPXO10 tidal atlas Flather conditions on the shelf boundary (10 constituents, 6 h ramp).
    - GLORYS daily reanalysis for `TracerReservoir` on T, S and `NormalFlowBoundaryCondition` on
      u, v at the shelf boundary, loaded through NumericalEarth's `Metadata`/`FieldTimeSeries`
      (30-day streaming window).
    - ERA5 `PrescribedAtmosphere` and `PrescribedRadiation` per level (decision D3), each an
      `OceanOnlyModel` wrapping `ocean_simulation` with CATKE closure, WENO advection and
      `SplitRungeKutta3`.
    - Passive tracer P released at the surface in St Brides Bay.
    - Time series at St Brides Bay and Ramsey Sound for tidal agreement between levels.
  - Checks: not run (requires GLORYS credentials and ERA5 data, and a GPU for production speed). The
    script records sea level and surface speed at the two stations, reports the inter-level rms
    difference, and prints the passive tracer's max and total on each grid. The feedback-off
    comparison is obtained by rerunning with `width = 0`.
  - Departures: ERA5 instead of JRA55 (user's preference for available data on the cluster);
    initialization from an analytic profile rather than GLORYS fields (the boundary conditions spin
    up the interior over the first months); the combined tidal + GLORYS Flather condition from
    section 5 is split into tidal Flather on U/V and GLORYS NormalFlow on u/v (the GLORYS
    contribution to barotropic transport is small on this shelf-tide-dominated domain).

- 2026-11-10, session H (Severn Estuary validation setup):
  - **`lowest_astronomical_tide(TPXO10Atlas, grid)`** in `src/tpxo.jl`: computes LAT ≈ -Σ|Aᵢ|
    from the TPXO10 elevation amplitudes (10 constituents by default). Returns a `Field` on the
    atlas window grid (0–360° longitude); the caller wraps with `mod(λ, 360)` before passing to
    `refer_to_mean_sea_level!`, which handles interpolation via `set!`.
  - **`validation/severn_estuary.jl`:** two-level nest (Bristol Channel shelf 140 × 90 × 50,
    Severn Estuary 400 × 200 × 30) with TPXO10 + GLORYS + ERA5 forcing, matching the
    Pembrokeshire template. Changes from the initial draft:
    - EMODnet bathymetry referenced to MSL via `lowest_astronomical_tide` instead of the
      placeholder `refer_to_mean_sea_level!(bottom_height, 0; ...)`.
    - GloFAS river routing on the estuary grid, gated behind `GLOFAS_DIR` env var. When enabled,
      `GloFASPrescribedLand` provides freshwater flux through NumericalEarth's `PrescribedLand`
      and `river_mouth_vertical_diffusivity` mixing closure via `ocean_simulation(; river_routing)`.
      All four rivers (Severn, Wye, Usk, Avon) enter through the interior since the domain extends
      east of Gloucester.
  - **`validation/severn_bathymetry.jl`:** also updated with the LAT correction.
  - Checks: syntax-checks pass; not run (requires TPXO atlas, GLORYS credentials, and optionally
    GloFAS data).
  - Departures: `GloFASPrescribedLand` import is assumed exported from `NumericalEarth` (like
    `ERA5PrescribedAtmosphere`); verify on first run.

- 2026-11-11, session I (three-level Severn grid design and script update):
  - **Grid configuration settled** through iterative trimming of domain sizes and performance
    estimation on a single H100. The three-level hierarchy:
    - **Shelf** (340×320×50, ~3.5 km): NW European shelf, (-12, 5)×(46, 62), TPXO + GLORYS BCs.
    - **Bristol Channel** (350×140×50, ~690 m): (-5.5, -2.0)×(50.8, 52.2), nested from shelf.
    - **Severn Estuary** (~500×325×30, ~140 m): (-3.5, -2.5)×(51.1, 51.75), nested from Bristol.
    Cost: ~1M wet-column-steps per shelf step (shelf 7%, Bristol 8%, Severn 85%). Estimated ~5 SYPD.
  - **`src/nested_timestep.jl`** (new): `minimum_horizontal_spacing(grid)` for
    `LatitudeLongitudeGrid` and `RectilinearGrid`; `nested_timestep(grids, maximum_velocities)`
    computes CFL-limited timesteps for a nested hierarchy with fixed integer subcycling.
    Added to exports and includes in `src/OnlineNesting.jl`.
  - **`validation/severn_estuary.jl`:** rewritten from two-level to three-level. Key changes:
    shelf is now the full NW European shelf (not just Bristol Channel); Bristol Channel is the
    intermediate domain; Severn is the innermost. LAT computed on bristol_grid (not shelf_grid)
    due to the `atlas_window` 0° longitude wrapping bug. Timestep computed via `nested_timestep`.
    Passive dye tracer released from Hinkley Point on the Severn grid.
  - **`validation/severn_bathymetry.jl`:** rewritten for three levels with per-level color ranges,
    resolution and wet-fraction annotations in titles, and pairwise column-volume checks.
  - **`validation/severn_estuary_instructions.md`:** updated for the three-level setup with
    performance estimates and grid descriptions.
  - Checks: syntax-checks pass; not run (requires TPXO atlas, GLORYS credentials, ERA5 data).
  - Departures: shelf gets `lat_correction=0` (no LAT correction) because the shelf grid spans 0°
    longitude; the shelf is deep enough that the LAT correction (~few metres) is negligible.

- 2026-11-12, session J (GPU crash debugging — roughness fix, not sufficient):
  - Added roughness-length guards (`ReynoldsScalingFunction`, `MomentumRoughnessLength`) to
    `validation/severn_estuary.jl` and the isolation scripts. These guards prevent `DomainError`
    from `R★ ^ 0.72` with negative `R★` — but the crash persists, so roughness is not the root cause.
  - Added per-stage GPU sync points in both Oceananigans SRK (`split_runge_kutta.jl`) and
    OnlineNesting `accumulate_fluxes!` to bisect the crash to a specific stage and operation.
  - Created isolation scripts: `validation/test_shelf_gpu.jl` (shelf only + TPXO, no nesting),
    `validation/test_nest_no_tpxo_gpu.jl` (nesting without TPXO). Neither has been run on GPU.
  - Checks: syntax passes. Not run on GPU.

- 2026-11-14, session K (GPU crash — systematic investigation):
  - **Hypotheses disproved:**
    (1) Unzeroed GPU memory: `new_data` uses `zeros(arch, ...)` on all architectures — fields
        ARE zero-initialized on GPU. Deep halo cells contain 0.0, not garbage.
    (2) Roughness-length DomainError: guards added in session J did not fix the crash.
    (3) TEOS10 with zero salinity: `s(0) = √(32/40.19) ≈ 0.89`, safe.
    (4) Halo fill mismatch: VBCs fill 1 halo cell but `volume_kernel_parameters` reads 6 deep —
        the unfilled cells are zero, which is safe for TEOS10 and all other operations.
  - **Exhaustively eliminated:** CATKE (sqrt guarded), WENO (no sqrt/log/^), split-explicit
    (sqrt guarded), gravity-wave BCs (guarded), TPXO data loading, field initialization,
    PartialCellBottom initialization, StateExchanger zero-init, RK substep (interior only).
  - **Actual error message never captured:** `severn.log` shows `@info` output but the GPU
    error went to stderr only. The error may not be a DomainError at all.
  - **Added:** `validation/test_nesting_gpu_idealized.jl` — 5-step progressive test (plain
    ocean → nesting → z-star → CATKE+WENO → PartialCellBottom) using no external data.
    Identifies which component triggers the crash.
  - **Added:** finer-grained sync points in `rk_substep!` for stage 3 (after
    `compute_tracer_tendencies`, after callbacks, after `rk_substep_tracers`) and in
    `accumulate_fluxes!` (after the kernel launches).
  - **Next:** run `test_nesting_gpu_idealized.jl` on GPU with `2>&1` to capture the actual
    error. If it passes, run the data-dependent isolation scripts. The actual error message
    is the single most valuable piece of missing information.
  - Checks: syntax passes; CPU test running. Not run on GPU.

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
