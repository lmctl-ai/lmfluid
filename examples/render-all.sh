#!/usr/bin/env bash
# Regenerate every docs/images/* file from the examples. Each example writes
# images/<name>.png (and images/<name>.gif for animated ones when ffmpeg is on PATH);
# this script runs them all against build/lmfluid and copies the results into docs/.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$PWD/build:$PATH"
compiler="${LMSCRIPT_BIN:-$(scripts/select-compiler.sh)}"
mkdir -p docs/images
command -v ffmpeg >/dev/null 2>&1 || echo "SKIP animations: ffmpeg unavailable (PNGs are still regenerated)"
for program in examples/*.lms; do
  name="$(basename "$program" .lms)"
  rm -f "images/$name.gif"
  "$compiler" "$program"
  cp "images/$name.png" "docs/images/$name.png"
  if [ -f "images/$name.gif" ]; then cp "images/$name.gif" "docs/images/$name.gif"; fi
done
python3 scripts/check-gallery.py --generated
