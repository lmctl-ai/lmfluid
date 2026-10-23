#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
BUILD_DIR="${BUILD_DIR:-build}"
HOST="$BUILD_DIR/lmfluid"
UNAME_S="$(uname -s)"
case "$UNAME_S" in MINGW*|MSYS*|CYGWIN*) HOST="$HOST.exe";; esac
WORK="$(mktemp -d "${TMPDIR:-/tmp}/lmfluid-test-XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
# Native Windows tools (node, the host) cannot address git-bash /tmp or /d/a paths.
case "$UNAME_S" in MINGW*|MSYS*|CYGWIN*) WORK="$(cd "$WORK" && pwd -W | tr '\134' '/')";; esac
pass=0
ok(){ pass=$((pass+1)); echo "ok $pass - $1"; }
fail(){ echo "FAIL: $*" >&2; exit 1; }
sha(){ if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -c1-64; else shasum -a 256 "$1" | cut -c1-64; fi; }

# lmscript resolves cli="lmfluid" on PATH first, so programs bind to the binary under test.
PATH="$(cd "$BUILD_DIR" 2>/dev/null && pwd || echo "$HERE/$BUILD_DIR"):$PATH"
export PATH
LMSCRIPT="$(scripts/select-compiler.sh)"
scripts/sync-sdk.sh --verify >/dev/null || fail "SDK header pin mismatch"
ok "SDK header matches sdk.lock"
scripts/sync-deps.sh --verify >/dev/null || fail "third_party pin mismatch"
ok "third_party matches deps.lock"

GEN=""
if [ ! -f "$BUILD_DIR/CMakeCache.txt" ] && command -v ninja >/dev/null 2>&1; then GEN="-G Ninja"; fi
# shellcheck disable=SC2086
cmake -S . -B "$BUILD_DIR" $GEN -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "$BUILD_DIR" --parallel "${BUILD_JOBS:-4}" >"$WORK/build.log" 2>&1 || { tail -40 "$WORK/build.log" >&2; fail "build"; }
ok "build"
ctest --test-dir "$BUILD_DIR" --output-on-failure || fail "native fluid tests"
ok "native fluid regression tests"
want_version="lmfluid $(sed -n 's/^project(lmfluid VERSION \([0-9.]*\).*/\1/p' CMakeLists.txt)"
[ "$("$HOST" --version)" = "$want_version" ] || fail "--version is not exactly '$want_version'"
ok "--version prints exactly '$want_version'"
python3 scripts/gen-schema.py >/dev/null && git diff --quiet -- schema/lmfluid.dsl.json 2>/dev/null || [ ! -e .git ] || fail "schema/lmfluid.dsl.json is stale: run scripts/gen-schema.py"
ok "schema matches the registrations in cpp/main.cpp"
"$HOST" --dsl-schema --json >"$WORK/schema.json" || fail "embedded schema print"
cmp -s "$WORK/schema.json" schema/lmfluid.dsl.json || fail "embedded schema differs from schema/lmfluid.dsl.json"
ok "--dsl-schema --json prints byte-identical embedded schema bytes"
diff -u <("$LMSCRIPT" --dsl-signatures schema/lmfluid.dsl.json) <("$HOST" --dsl-signatures) || fail "host/compiler signatures differ"
ok "host registry matches compiler signatures"
"$HOST" --notices >"$WORK/notices.txt" || fail "--notices"
grep -q "MFEM" "$WORK/notices.txt" || fail "--notices missing third-party text"
cmp -s "$WORK/notices.txt" THIRD_PARTY_NOTICES || fail "--notices differs from THIRD_PARTY_NOTICES"
ok "--notices prints THIRD_PARTY_NOTICES byte for byte"

case "$UNAME_S" in
  Linux) if ldd "$HOST" >/dev/null 2>&1; then fail "executable is dynamically linked"; fi; ok "executable is fully static (no dynamic loader)";;
  Darwin) bad="$(otool -L "$HOST" | tail -n +2 | grep -v -E '/usr/lib/(libSystem|libc\+\+|libresolv)' || true)"
          [ -z "$bad" ] || fail "non-system dylib dependencies: $bad"; ok "only system libraries linked";;
  MINGW*|MSYS*|CYGWIN*)
    if command -v objdump >/dev/null 2>&1; then
      bad="$(objdump -p "$HOST" | grep 'DLL Name' | grep -i -E 'libstdc|libgcc|winpthread|libgomp' || true)"
      [ -z "$bad" ] || fail "MinGW runtime DLLs are imported: $bad"
      ok "no MinGW runtime DLL imports"
    else ok "static-link check skipped (no objdump)"; fi;;
  *) ok "static-link check skipped on $UNAME_S";;
esac

# Programs: compile every test and example, run them (they assert their own numbers).
mkdir -p build/examples images
rm -f images/*.png
n=0
for src in tests/*.lms examples/*.lms; do
  name="$(basename "$src" .lms)"
  [ ! -e "$WORK/compiled-$name" ] || fail "duplicate program basename: $name"
  touch "$WORK/compiled-$name"
  "$LMSCRIPT" "$src" "build/$name.dslb" >"$WORK/$name.log" 2>&1 || { cat "$WORK/$name.log" >&2; fail "$name compile"; }
  n=$((n+1))
done
ok "$n programs compile as bytecode v5"
for src in tests/*.lms examples/*.lms; do
  name="$(basename "$src" .lms)"
  "$HOST" --dsl "build/$name.dslb" >"$WORK/$name.out" 2>"$WORK/$name.err" || { cat "$WORK/$name.out" "$WORK/$name.err" >&2; fail "$name run"; }
  grep -q '"ok": true' "$WORK/$name.out" || { cat "$WORK/$name.out" >&2; fail "$name did not report ok"; }
done
ok "$n programs run and pass their own assertions"

# Determinism: the CSV files the examples write are identical on rerun.
for f in build/examples/*.csv; do echo "$(sha "$f") $f"; done >"$WORK/before.txt"
for src in examples/*.lms; do
  name="$(basename "$src" .lms)"
  "$HOST" --dsl "build/$name.dslb" >/dev/null 2>&1 || fail "$name rerun"
done
for f in build/examples/*.csv; do echo "$(sha "$f") $f"; done >"$WORK/after.txt"
cmp -s "$WORK/before.txt" "$WORK/after.txt" || fail "result files differ on rerun"
ok "all $(wc -l <"$WORK/after.txt" | tr -d ' ') result files are byte-identical on rerun"

# Error paths: each snippet must fail the run with exit 3 and a message naming the problem.
expect_error() { # name, snippet, expected message fragment
  cat >"$WORK/$1.lms" <<LMS
p = _dsl_ cli="lmfluid";
$2
end _dsl_
savedsl(p, "$WORK/$1.dslb");
LMS
  "$LMSCRIPT" "$WORK/$1.lms" >/dev/null 2>"$WORK/$1.cerr" || { cat "$WORK/$1.cerr" >&2; fail "$1 compile"; }
  set +e; out="$("$HOST" --dsl "$WORK/$1.dslb" 2>/dev/null)"; code=$?; set -e
  [ "$code" = 3 ] || fail "$1: exit code $code, expected 3 ($out)"
  echo "$out" | grep -q "$3" || fail "$1: message lacks '$3': $out"
}
expect_error grid 'let m = fluid.rect(0,4,1.0,1.0);' "mesh"
expect_error dimension 'let m = fluid.rect(4,4,-1.0,1.0);' "mesh"
expect_error order 'let p = fluid.poisson(fluid.rect(2,2,1.0,1.0),0);' "order"
expect_error diffusion 'let f = fluid.poisson(fluid.rect(4,4,1.0,1.0),1).diffusion("-1").solve();' "diffusion"
expect_error expression 'let f = fluid.poisson(fluid.rect(4,4,1.0,1.0),1).source("unknown+1").solve();' "unknown"
expect_error outside 'let f = fluid.poisson(fluid.rect(4,4,1.0,1.0),1).solve(); let x = f.value(2.0,2.0);' "outside"
expect_error missing_mesh 'let m = fluid.load_mesh("missing-fluid.mesh");' "mesh"
expect_error heat_step 'let f = fluid.heat(fluid.rect(4,4,1.0,1.0),1).evolve(1.0,0.0);' "step"
expect_error viscosity 'let f = fluid.channel(8,8,0.0,0.00001);' "viscosity"
expect_error lid 'let f = fluid.cavity(8,8,0.1,0.5);' "lid"
expect_error velocity_index 'let f = fluid.periodic(8,8,0.1); let x=f.velocity(100,0,0);' "index"
expect_error matrix_index 'let m=fluid.periodic(8,8,0.1).sample(); let x=m.get(10000,0);' "index"
expect_error assertfail 'fluid.assert_close(1.0,2.0,0.1);' "assert_close"
expect_error nan_assert 'fluid.assert_close(sqrt(-1.0),0.0,0.1);' "assert_close"
expect_error assert_true 'fluid.assert_true(false);' "assert_true"
# Malformed MFEM meshes are rejected before MFEM's reader (which segfaulted on a dangling vertex index).
printf 'MFEM mesh v1.0\n\ndimension\n2\n\nelements\n1\n1 3 0 1 2 3\n\nboundary\n0\n\nvertices\n2\n2\n0 0\n1 0\n' >"$WORK/dangling.mesh"
expect_error dangling_vertex "let m = fluid.load_mesh(\"$WORK/dangling.mesh\");" "only 2 vertices"
printf 'garbage\n' >"$WORK/garbage.mesh"
expect_error mesh_format "let m = fluid.load_mesh(\"$WORK/garbage.mesh\");" "MFEM mesh v1.0"
expect_error no_convergence 'let f = fluid.poisson(fluid.rect(16,16,1.0,1.0),1).source("1").tol(0.000000000001, 1).solve();' "did not converge"
[ "$(printf '%s\n' "$out" | wc -l | tr -d ' ')" = 1 ] || fail "solver diagnostics leaked into the JSON stdout: $out"
# Expression strings: deep nesting used to overflow the parser's stack (SIGSEGV).
deep="$(printf '(%.0s' $(seq 300))1$(printf ')%.0s' $(seq 300))"
expect_error deep_expression "let x = fluid.eval(\"$deep\");" "nested too deeply"
expect_error huge_number 'let x = fluid.eval("1e400");' "out of range"
ok "errors stop with exit 3 and clear messages (20 cases)"

# Dry run loads/verifies bytecode without executing numerical work or writing CSV.
rm -f build/examples/poisson.csv
"$HOST" --dry-run --dsl build/poisson.dslb >"$WORK/dry.out" || fail "dry run"
[ ! -e build/examples/poisson.csv ] || fail "dry run created CSV"
"$HOST" --dsl build/poisson.dslb >/dev/null || fail "restore example CSV"
ok "dry run verifies without writing outputs"

for flags in '--unknown' '--version extra' '--dsl-schema --bogus' '--dry-run' '--dsl --bogus'; do
  set +e
  # Deliberate word splitting: each fixed test string is a CLI argument sequence.
  # shellcheck disable=SC2086
  "$HOST" $flags >/dev/null 2>&1
  code=$?
  set -e
  [ "$code" = 2 ] || fail "usage flags '$flags': exit $code"
done
ok "invalid CLI combinations use exit 2"

set +e; "$HOST" --dsl "$WORK/does-not-exist.dslb" >/dev/null 2>&1; code=$?; set -e
[ "$code" = 3 ] || fail "missing bytecode exit code $code, expected 3"
set +e; "$HOST" >/dev/null 2>&1; code=$?; set -e
[ "$code" = 2 ] || fail "usage exit code $code, expected 2"
ok "usage and load errors use exit codes 2 and 3"

python3 scripts/check-gallery.py --generated
ok "gallery links, image signatures, dimensions and byte budgets"

echo "all $pass checks passed"
