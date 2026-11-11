# Running the Severn Estuary three-level nesting case on bouchet

## 1. Clone and set up Oceananigans

The simulation requires the open-boundary filtered-transport fix (PR #6485), which is not in
released Oceananigans 0.117.2. Check out the PR branch and point `Project.toml` at it:

```bash
cd ~/Documents
git clone https://github.com/CliMA/Oceananigans.jl.git Oceananigans-open-boundary-fix
cd Oceananigans-open-boundary-fix
git checkout jsw/open-boundary-filtered-transport
```

Then edit `OnlineNesting.jl/Project.toml` so the `[sources]` path points to the checkout:

```toml
[sources]
Oceananigans = {path = "/home/jago/Documents/Oceananigans-open-boundary-fix"}
```

If the branch fails to resolve (e.g. a DocStringExtensions compat mismatch), fall back to
cherry-picking the fix onto the release tag:

```bash
git checkout v0.117.2
git cherry-pick <commit-hash-from-PR-6485>
```

## 2. Data

| Dataset | How to get it | Env var |
|---------|--------------|---------|
| TPXO10 atlas | Copy from laptop or request at https://www.tpxo.net | `TPXO_DIR` (default `~/Downloads/TPXO10_atlas_v2_nc/`) |
| GLORYS daily | Download script or auto-downloads on first run | `NUMERICALEARTH_DATA_DIRECTORY` |
| ERA5 hourly | Auto-downloads on first run | `ERA5_DIR` (point to a shared cache) |
| GloFAS reanalysis | Download script (optional; omit to skip rivers) | `NUMERICALEARTH_DATA_DIRECTORY` |

### Download GLORYS and GloFAS

The download script fetches GLORYS reanalysis (T, S, u, v, η for the shelf boundary) and
GloFAS river discharge (Severn, Wye, Usk, Avon) via NumericalEarth's data extensions. Both
`CopernicusMarine` and `CDSAPI` are Julia packages — no system tools needed.

Set up credentials first:

```bash
# Copernicus Marine (GLORYS) — run once:
julia --project=test -e 'using CopernicusMarine; CopernicusMarine.login()'

# CDS API (GloFAS) — create ~/.cdsapirc with your key from
# https://cds.climate.copernicus.eu/how-to-api
```

Then download (set `NUMERICALEARTH_DATA_DIRECTORY` to a shared path on the cluster):

```bash
NUMERICALEARTH_DATA_DIRECTORY=/path/to/data \
    julia --project=test validation/download_severn_data.jl
```

Already-downloaded files are skipped, so it can be rerun after interruptions.

### TPXO10

Copy the atlas directory from your laptop. The required files are `grid_tpxo10atlas_v2.nc` plus
`h_<constituent>_tpxo10_atlas_30_v2.nc` and `u_<constituent>_tpxo10_atlas_30_v2.nc` for each of
M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm:

```bash
rsync -av ~/Downloads/TPXO10_atlas_v2_nc/ bouchet:Downloads/TPXO10_atlas_v2_nc/
```

### ERA5

ERA5 is downloaded automatically when the simulation runs. Point to a shared cache if available:

```bash
export ERA5_DIR="/path/to/era5/cache"
```

## 3. Environment

Instantiate the test environment (which has CopernicusMarine, CDSAPI and ERA5 extensions):

```bash
cd /path/to/OnlineNesting.jl
julia --project=test -e 'using Pkg; Pkg.instantiate()'
```

## 4. Run

### CPU (for testing the setup — start with a short run)

```bash
julia --project=test -t auto validation/severn_estuary.jl 7   # 7 days
```

### GPU (for production)

Change `arch = CPU()` to `arch = GPU()` at line 33 of `severn_estuary.jl`, then:

```bash
julia --project=test -t auto -e 'using CUDA; include("validation/severn_estuary.jl")'
```

The first positional argument is the number of days to run (default 30). Start with 7 to check
everything loads, then run 30+ for tidal analysis.

## 5. What it does

- Builds three nested grids with 5× refinement between levels:
  - **Shelf** (340×320×50, ~3.5 km): NW European shelf with TPXO tidal + GLORYS open boundaries
  - **Bristol Channel** (350×140×50, ~690 m): nested from shelf, covers the Bristol Channel
  - **Severn Estuary** (~500×325×30, ~140 m): nested from Bristol, Llantwit Major to Gloucester
- Computes timesteps from CFL constraints using `nested_timestep` (RK3, CFL ≈ 1.7)
- Downloads EMODnet bathymetry, corrects from LAT to MSL using TPXO10 tidal amplitudes (computed
  on the Bristol grid to avoid the 0° longitude atlas_window bug), blends between levels
- Creates an ERA5-forced `EarthSystemModel` per level with CATKE closure
- Optionally routes GloFAS river discharge onto the Severn grid via `PrescribedLand`
- Records sea level and surface speed at Avonmouth, Newport, Mumbles, Ilfracombe and Hinkley Point
- Reports inter-level agreement (shelf↔Bristol, Bristol↔Severn) at the end

## 6. Performance

Estimated ~5 SYPD on a single H100. The Severn domain dominates cost (~85% of wet-column-steps).
A 30-day run should complete well within an hour on GPU.

## 7. What to look for

- **Tidal range agreement**: adjacent levels should produce similar tidal ranges where they
  overlap. The Bristol Channel has some of the largest tides in the world (~14 m spring range at
  Avonmouth), so even small inconsistencies are visible.
- **Severn salinity** (if GloFAS is enabled): river freshwater should create a salinity gradient
  from ~0 psu at the tidal limit to ~35 psu in the outer Bristol Channel.
- **Spin-up**: the boundary conditions spin up the interior over the first few days. Initial
  conditions are an analytic profile, not GLORYS fields.

## 8. Feedback-off comparison

To compare with one-way nesting, change the `TwoWayNestedModel` calls near the bottom of the
script to use `width = 0`.
