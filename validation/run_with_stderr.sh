#!/bin/bash
# Run a Julia script with stderr interleaved into stdout, so the error message
# appears in the log alongside the @info output.
#
# Usage: bash validation/run_with_stderr.sh validation/severn_estuary.jl 0.1
#    or: bash validation/run_with_stderr.sh validation/test_nesting_gpu_idealized.jl

SCRIPT="${1:?Usage: $0 <script.jl> [args...]}"
shift

julia --project=test -t auto "$SCRIPT" "$@" 2>&1 | tee "${SCRIPT%.jl}.combined.log"
