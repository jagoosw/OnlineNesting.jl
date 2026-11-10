# Running the Pembrokeshire three-level nesting case

## 1. Data

Three datasets are needed:

| Dataset | How to get it | Where it lives |
|---------|--------------|----------------|
| TPXO10 atlas | Already available | `~/Downloads/TPXO10_atlas_v2_nc/` |
| GLORYS daily | Download via Copernicus Marine | Auto-cached by NumericalEarth |
| ERA5 hourly | Already on bouchet, or auto-downloads | Set `ERA5_DIR` to the shared cache |

For GLORYS, set the Copernicus Marine credentials before running:

```bash
export COPERNICUSMARINE_SERVICE_USERNAME="your-username"
export COPERNICUSMARINE_SERVICE_PASSWORD="your-password"
```

NumericalEarth downloads GLORYS and ERA5 on first use through its `FieldTimeSeries`/`Metadata`
system. If ERA5 is already on bouchet, point to it:

```bash
export ERA5_DIR="/path/to/era5/cache"
```

## 2. Environment

The script needs NumericalEarth's test environment (for the CopernicusMarine extension and
ERA5/GLORYS dataset types):

```bash
julia --project=test -e 'using Pkg; Pkg.instantiate()'
```

## 3. Run

On CPU (slow, for testing the setup):

```bash
julia --project=test -t auto validation/pembrokeshire.jl 7   # 7 days
```

On GPU (for production), change `arch = CPU()` to `arch = GPU()` at line 28 of the script, then:

```bash
julia --project=test -t auto -e 'using CUDA; include("validation/pembrokeshire.jl")'
```

The first argument is the number of days (default 730, i.e. 2 years). Start with a short run to
check everything loads.

## 4. What it does

- Builds shelf (220×245), approaches (250×245) and bay (305×290) grids with 50 z-star layers
- Loads EMODnet bathymetry, reconciles between levels
- Applies TPXO10 tidal Flather conditions and GLORYS T/S/u/v at the shelf boundary
- Creates an ERA5-forced `EarthSystemModel` per level with CATKE closure
- Releases a passive tracer at the surface in St Brides Bay
- Records sea level and surface speed at St Brides Bay and Ramsey Sound
- Reports inter-level agreement at the end

## 5. Feedback-off comparison

To compare with one-way nesting, change the `TwoWayNestedModel` calls near the bottom of the
script to use `width = 0`.
