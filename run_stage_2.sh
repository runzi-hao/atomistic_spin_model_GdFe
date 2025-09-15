#!/usr/bin/env bash
set -euo pipefail

COMMIT="${GIT_COMMIT:?Must pass GIT_COMMIT env var}"
OUTROOT="${OUTROOT:-outputs}"
BUILD_CMD="${BUILD_CMD:-make -j}"
RUN_CMD="${RUN_CMD:-./my_program}"

# 1) Sync refs, ensure commit exists on origin
git remote update --prune
if ! git cat-file -e "$COMMIT^{commit}" 2>/dev/null; then
  echo "Commit $COMMIT not found locally. Fetching full refs..."
  git fetch --all --prune
  git cat-file -e "$COMMIT^{commit}"  # will fail if truly missing
fi

# 2) Checkout EXACT commit (detached HEAD)
git checkout -q --detach "$COMMIT"
[[ "$(git rev-parse HEAD)" == "$COMMIT" ]]

# 3) Stamp outputs
STAMP="$(date +%Y%m%d_%H%M%S)"
DESC="$(git describe --always --dirty --tags 2>/dev/null || git rev-parse --short HEAD)"
OUTDIR="${OUTROOT}/${DESC}_${STAMP}"
mkdir -p "$OUTDIR"

{
  echo "commit: $(git rev-parse HEAD)"
  echo "describe: $DESC"
  echo "remote: origin ($(git remote get-url origin))"
  echo "date: $(date --iso-8601=seconds)"
} > "$OUTDIR/version.txt"

# 4) Build & run; capture logs
( $BUILD_CMD ) > "$OUTDIR/build.log" 2>&1
( $RUN_CMD   ) > "$OUTDIR/run.log"   2>&1

echo "DONE → $OUTDIR"
