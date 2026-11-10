# Running the Severn Estuary two-level nesting case on bouchet

## 1. Clone and set up Oceananigans

The simulation requires the open-boundary filtered-transport fix (PR #6485), which is not in
released Oceananigans 0.117.2. On your laptop the `[sources]` pin in `Project.toml` points to a
local patched copy; on bouchet you need to recreate it:

```bash
cd ~/Documents
git clone https://github.com/CliMA/Oceananigans.jl.git Oceananigans-0.117.2-open-boundary-fix
cd Oceananigans-0.117.2-open-boundary-fix
git checkout v0.117.2
# Apply the filtered-transport fix from PR #6485:
git cherry-pick <commit-hash>   # or apply the diff manually
```

Then edit `OnlineNesting.jl/Project.toml` so the `[sources]` path points to wherever you put it
on bouchet:

```toml
[sources]
Oceananigans = {path = "/home/jago/Documents/Oceananigans-0.117.2-open-boundary-fix"}
```

## 2. Data

| Dataset | How to get it | Env var |
|---------|--------------|---------|
| TPXO10 atlas | Copy from laptop or request at https://www.tpxo.net | `TPXO_DIR` (default `~/Downloads/TPXO10_atlas_v2_nc/`) |
| GLORYS daily | Auto-downloads via Copernicus Marine credentials | `GLORYS_DIR` (optional cache dir) |
| ERA5 hourly | Already on bouchet, or auto-downloads | `ERA5_DIR` (point to the shared cache) |
| GloFAS reanalysis | Optional; download via CDS API or copy | `GLOFAS_DIR` (omit to skip rivers) |

### TPXO10

Copy the atlas directory from your laptop. The required files are `grid_tpxo10atlas_v2.nc` plus
`h_<constituent>_tpxo10_atlas_30_v2.nc` and `u_<constituent>_tpxo10_atlas_30_v2.nc` for each of
M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm:

```bash
rsync -av ~/Downloads/TPXO10_atlas_v2_nc/ bouchet:Downloads/TPXO10_atlas_v2_nc/
```

### GLORYS

Set Copernicus Marine credentials before running:

```bash
export COPERNICUSMARINE_SERVICE_USERNAME="your-username"
export COPERNICUSMARINE_SERVICE_PASSWORD="your-password"
```

NumericalEarth downloads GLORYS on first use. To cache it:

```bash
export GLORYS_DIR="/path/to/glorys/cache"
```

### ERA5

If ERA5 is already on bouchet:

```bash
export ERA5_DIR="/path/to/era5/cache"
```

### GloFAS (optional — river freshwater forcing)

GloFAS provides freshwater discharge for the Severn, Wye, Usk and Avon rivers. Without it the
simulation still runs but has no river forcing. To enable it, download GloFAS reanalysis data
(needs EWDS credentials from https://cds.climate.copernicus.eu) and point to the cache:

```bash
export GLOFAS_DIR="/path/to/glofas/cache"
```

## 3. Environment

Instantiate the test environment (which has the CopernicusMarine and ERA5 extensions):

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

Change `arch = CPU()` to `arch = GPU()` at line 31 of `severn_estuary.jl`, then:

```bash
julia --project=test -t auto -e 'using CUDA; include("validation/severn_estuary.jl")'
```

The first positional argument is the number of days to run (default 30). Start with 7 to check
everything loads, then run 30+ for tidal analysis.

## 5. What it does

- Builds a shelf grid (140×90×50) and estuary grid (400×200×30, 5× refinement) with z-star layers
- Downloads EMODnet bathymetry, corrects from LAT to MSL using TPXO10 tidal amplitudes, blends
  between levels
- Applies TPXO10 tidal Flather conditions and GLORYS T/S/u/v at the shelf boundary
- Creates an ERA5-forced `EarthSystemModel` per level with CATKE closure
- Optionally routes GloFAS river discharge onto the estuary grid via `PrescribedLand`
- Records sea level and surface speed at Avonmouth, Newport, Mumbles and Ilfracombe
- Reports inter-level agreement at the end

## 6. What to look for

- **Tidal range agreement**: the shelf and estuary should produce similar tidal ranges where they
  overlap. The Bristol Channel has some of the largest tides in the world (~14 m spring range at
  Avonmouth), so even small inconsistencies are visible.
- **Estuary salinity** (if GloFAS is enabled): river freshwater should create a salinity gradient
  from ~0 psu at the tidal limit to ~35 psu in the outer Bristol Channel.
- **Spin-up**: the boundary conditions spin up the interior over the first few days. Initial
  conditions are an analytic profile, not GLORYS fields.

## 7. Feedback-off comparison

To compare with one-way nesting, change the `TwoWayNestedModel` call near the bottom of the
script to use `width = 0`.
