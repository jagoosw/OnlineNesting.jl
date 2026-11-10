#!/bin/bash
# Download GLORYS reanalysis data for the Pembrokeshire three-level nesting case.
#
# Prerequisites:
#   pip install copernicusmarine
#   export COPERNICUSMARINE_SERVICE_USERNAME="..."
#   export COPERNICUSMARINE_SERVICE_PASSWORD="..."
#
# TPXO10-atlas-v2 is licensed and cannot be downloaded here; request it at https://www.tpxo.net
# and unpack into $TPXO_DIR.
#
# ERA5 is downloaded automatically by NumericalEarth when the script runs.

set -euo pipefail

GLORYS_DIR="${GLORYS_DIR:-data/glorys}"
TPXO_DIR="${TPXO_DIR:-$HOME/Downloads/TPXO10_atlas_v2_nc}"
START_DATE="${START_DATE:-2020-01-01}"
END_DATE="${END_DATE:-2022-01-01}"

# Shelf domain with a 0.5° margin for interpolation
WEST=-11.0
EAST=-2.0
SOUTH=48.5
NORTH=55.0

mkdir -p "$GLORYS_DIR"

echo "Downloading GLORYS daily reanalysis for ${START_DATE} to ${END_DATE}..."
copernicusmarine subset \
    --dataset-id cmems_mod_glo_phy_my_0.083deg_P1D-m \
    --variable thetao --variable so --variable uo --variable vo --variable zos \
    --minimum-longitude "$WEST" --maximum-longitude "$EAST" \
    --minimum-latitude "$SOUTH" --maximum-latitude "$NORTH" \
    --start-datetime "${START_DATE}T00:00:00" \
    --end-datetime "${END_DATE}T00:00:00" \
    --output-directory "$GLORYS_DIR" \
    --output-filename glorys_pembrokeshire.nc \
    --force-download

echo ""
echo "GLORYS data saved to ${GLORYS_DIR}/glorys_pembrokeshire.nc"

echo ""
echo "Checking TPXO10 atlas files in ${TPXO_DIR}..."
MISSING=0
for FILE in grid_tpxo10atlas_v2.nc \
            h_m2_tpxo10_atlas_30_v2.nc u_m2_tpxo10_atlas_30_v2.nc \
            h_s2_tpxo10_atlas_30_v2.nc u_s2_tpxo10_atlas_30_v2.nc \
            h_n2_tpxo10_atlas_30_v2.nc u_n2_tpxo10_atlas_30_v2.nc \
            h_k2_tpxo10_atlas_30_v2.nc u_k2_tpxo10_atlas_30_v2.nc \
            h_k1_tpxo10_atlas_30_v2.nc u_k1_tpxo10_atlas_30_v2.nc \
            h_o1_tpxo10_atlas_30_v2.nc u_o1_tpxo10_atlas_30_v2.nc \
            h_p1_tpxo10_atlas_30_v2.nc u_p1_tpxo10_atlas_30_v2.nc \
            h_q1_tpxo10_atlas_30_v2.nc u_q1_tpxo10_atlas_30_v2.nc \
            h_mf_tpxo10_atlas_30_v2.nc u_mf_tpxo10_atlas_30_v2.nc \
            h_mm_tpxo10_atlas_30_v2.nc u_mm_tpxo10_atlas_30_v2.nc; do
    if [ ! -f "${TPXO_DIR}/${FILE}" ]; then
        echo "  MISSING: ${FILE}"
        MISSING=$((MISSING + 1))
    fi
done

if [ "$MISSING" -gt 0 ]; then
    echo "${MISSING} TPXO file(s) missing. Request access at https://www.tpxo.net"
    exit 1
else
    echo "  All TPXO constituent files present."
fi

echo ""
echo "ERA5 atmosphere and radiation data will be downloaded automatically by NumericalEarth."
echo "Set ERA5_DIR to a shared cache if available (e.g. on bouchet)."
echo ""
echo "Done. Run the case with:"
echo "  julia --project -t auto validation/pembrokeshire.jl"
