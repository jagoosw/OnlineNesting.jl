# Idealized GPU nesting test: no data dependencies, tiny grids, runs on any GPU.
# Exercises every OnlineNesting kernel: _restrict!, _accumulate_flux!, _correct_tracer!, _correct_sea_level!.
#
# Run: julia --project=test -t auto validation/test_nesting_gpu_idealized.jl
#
# Steps through progressively more complex setups. The first failure isolates the component.

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: Open, MutableVerticalDiscretization
using OnlineNesting

arch = GPU()

stratification(λ, φ, z) = 10 + z / 20

# ── Step 1: plain HydrostaticFreeSurfaceModel, no nesting ────────────────────

@info "Step 1: plain ocean model on GPU, one time step..."
grid1 = LatitudeLongitudeGrid(arch; size=(40, 40, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0),
                              halo=(7, 7, 7), topology=(Open, Open, Bounded))
grid1 = ImmersedBoundaryGrid(grid1, GridFittedBottom(-100))
model1 = HydrostaticFreeSurfaceModel(grid1; free_surface=SplitExplicitFreeSurface(grid1; substeps=20),
                                     timestepper=:SplitRungeKutta3, tracers=(:T, :S), buoyancy=SeawaterBuoyancy())
set!(model1, T=stratification, S=35)
time_step!(model1, 300)
@info "Step 1 passed."

# ── Step 2: two-level nest, GridFittedBottom, fixed z ────────────────────────

@info "Step 2: two-level nest on GPU, GridFittedBottom, fixed z..."
parent_grid = LatitudeLongitudeGrid(arch; size=(40, 40, 4), longitude=(0, 10), latitude=(40, 50), z=(-100, 0),
                                    halo=(7, 7, 7), topology=(Open, Open, Bounded))
parent_grid = ImmersedBoundaryGrid(parent_grid, GridFittedBottom(-100))
child_grid = nested_grid(parent_grid; longitude=(3, 7), latitude=(43, 47), refinement=3)

parent_model = HydrostaticFreeSurfaceModel(parent_grid; free_surface=SplitExplicitFreeSurface(parent_grid; substeps=20),
                                           timestepper=:SplitRungeKutta3, tracers=(:T, :S), buoyancy=SeawaterBuoyancy())
exchanger = StateExchanger(parent_model)
child_model = HydrostaticFreeSurfaceModel(child_grid;
                  free_surface=SplitExplicitFreeSurface(child_grid; substeps=20),
                  timestepper=:SplitRungeKutta3, tracers=(:T, :S), buoyancy=SeawaterBuoyancy(),
                  boundary_conditions=nested_boundary_conditions(exchanger, child_grid),
                  forcing=nested_forcings(exchanger, child_grid; rate=1/1hour, width=4))

nest = TwoWayNestedModel(parent_model, child_model, exchanger; width=2)
set!(parent_model, T=stratification, S=35)
set!(child_model, T=stratification, S=35)

simulation = Simulation(nest; Δt=300, stop_time=300)
run!(simulation)
@info "Step 2 passed."

# ── Step 3: z-star coordinate ────────────────────────────────────────────────

@info "Step 3: two-level nest on GPU, z-star (MutableVerticalDiscretization)..."
z_star = MutableVerticalDiscretization((-100, 0))

parent_grid3 = LatitudeLongitudeGrid(arch; size=(40, 40, 4), longitude=(0, 10), latitude=(40, 50), z=z_star,
                                     halo=(7, 7, 7), topology=(Open, Open, Bounded))
parent_grid3 = ImmersedBoundaryGrid(parent_grid3, GridFittedBottom(-100))
child_grid3 = nested_grid(parent_grid3; longitude=(3, 7), latitude=(43, 47), refinement=3)

parent_model3 = HydrostaticFreeSurfaceModel(parent_grid3; free_surface=SplitExplicitFreeSurface(parent_grid3; substeps=20),
                                            timestepper=:SplitRungeKutta3, tracers=(:T, :S), buoyancy=SeawaterBuoyancy())
exchanger3 = StateExchanger(parent_model3)
child_model3 = HydrostaticFreeSurfaceModel(child_grid3;
                   free_surface=SplitExplicitFreeSurface(child_grid3; substeps=20),
                   timestepper=:SplitRungeKutta3, tracers=(:T, :S), buoyancy=SeawaterBuoyancy(),
                   boundary_conditions=nested_boundary_conditions(exchanger3, child_grid3),
                   forcing=nested_forcings(exchanger3, child_grid3; rate=1/1hour, width=4))

nest3 = TwoWayNestedModel(parent_model3, child_model3, exchanger3; width=2)
set!(parent_model3, T=stratification, S=35)
set!(child_model3, T=stratification, S=35)

simulation3 = Simulation(nest3; Δt=300, stop_time=300)
run!(simulation3)
@info "Step 3 passed."

# ── Step 4: CATKE closure + WENO(order=7) + z-star ───────────────────────────

@info "Step 4: two-level nest on GPU, CATKE + WENO7 + z-star..."
using Oceananigans.Advection: WENO
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

vi = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)
weno = WENO(order=7, time_discretization=vi, minimum_buffer_upwind_order=1)

parent_grid4 = LatitudeLongitudeGrid(arch; size=(40, 40, 10), longitude=(0, 10), latitude=(40, 50),
                                     z=MutableVerticalDiscretization((-100, 0)),
                                     halo=(7, 7, 7), topology=(Open, Open, Bounded))
parent_grid4 = ImmersedBoundaryGrid(parent_grid4, GridFittedBottom(-100))
child_grid4 = nested_grid(parent_grid4; longitude=(3, 7), latitude=(43, 47), refinement=3)

parent_model4 = HydrostaticFreeSurfaceModel(parent_grid4;
                    free_surface=SplitExplicitFreeSurface(parent_grid4; substeps=20),
                    timestepper=:SplitRungeKutta3, tracers=(:T, :S), buoyancy=SeawaterBuoyancy(),
                    tracer_advection=weno, closure=CATKEVerticalDiffusivity())
exchanger4 = StateExchanger(parent_model4)
child_model4 = HydrostaticFreeSurfaceModel(child_grid4;
                   free_surface=SplitExplicitFreeSurface(child_grid4; substeps=20),
                   timestepper=:SplitRungeKutta3, tracers=(:T, :S), buoyancy=SeawaterBuoyancy(),
                   tracer_advection=weno, closure=CATKEVerticalDiffusivity(),
                   boundary_conditions=nested_boundary_conditions(exchanger4, child_grid4),
                   forcing=nested_forcings(exchanger4, child_grid4; rate=1/1hour, width=4))

nest4 = TwoWayNestedModel(parent_model4, child_model4, exchanger4; width=2)
set!(parent_model4, T=stratification, S=35)
set!(child_model4, T=stratification, S=35)

simulation4 = Simulation(nest4; Δt=300, stop_time=300)
run!(simulation4)
@info "Step 4 passed."

# ── Step 5: PartialCellBottom with realistic coast ───────────────────────────

@info "Step 5: two-level nest on GPU, PartialCellBottom + coast..."
coast(λ, φ) = ifelse(λ < 3, 0, ifelse(φ < 43.5, -50, -100))

parent_grid5 = LatitudeLongitudeGrid(arch; size=(40, 40, 10), longitude=(0, 10), latitude=(40, 50),
                                     z=MutableVerticalDiscretization((-100, 0)),
                                     halo=(7, 7, 7), topology=(Open, Open, Bounded))
parent_grid5 = ImmersedBoundaryGrid(parent_grid5, PartialCellBottom(coast))
child_grid5 = nested_grid(parent_grid5; longitude=(3, 7), latitude=(43, 47), refinement=3)

parent_model5 = HydrostaticFreeSurfaceModel(parent_grid5;
                    free_surface=SplitExplicitFreeSurface(parent_grid5; substeps=20),
                    timestepper=:SplitRungeKutta3, tracers=(:T, :S), buoyancy=SeawaterBuoyancy(),
                    tracer_advection=weno, closure=CATKEVerticalDiffusivity())
exchanger5 = StateExchanger(parent_model5)
child_model5 = HydrostaticFreeSurfaceModel(child_grid5;
                   free_surface=SplitExplicitFreeSurface(child_grid5; substeps=20),
                   timestepper=:SplitRungeKutta3, tracers=(:T, :S), buoyancy=SeawaterBuoyancy(),
                   tracer_advection=weno, closure=CATKEVerticalDiffusivity(),
                   boundary_conditions=nested_boundary_conditions(exchanger5, child_grid5),
                   forcing=nested_forcings(exchanger5, child_grid5; rate=1/1hour, width=4))

nest5 = TwoWayNestedModel(parent_model5, child_model5, exchanger5; width=2)
set!(parent_model5, T=stratification, S=35)
set!(child_model5, T=stratification, S=35)

simulation5 = Simulation(nest5; Δt=300, stop_time=300)
run!(simulation5)
@info "Step 5 passed."

@info "All steps passed — nesting machinery works on GPU with this configuration."
