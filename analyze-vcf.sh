#!/usr/bin/env bash
# analyze-vcf.sh — compatibility shim. The analyzer now lives in analyze.sh,
# which auto-detects snp-indel / sv / cnv. This forwards all arguments.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/analyze.sh" "$@"
