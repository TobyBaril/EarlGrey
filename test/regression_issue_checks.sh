#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

assert_contains() {
    local file="$1"
    local pattern="$2"
    local message="$3"
    if ! grep -Eq "$pattern" "$file"; then
        echo "FAIL: $message"
        exit 1
    fi
}

for script in earlGrey earlGreyLibConstruct earlGreyAnnotationOnly; do
    assert_contains "$script" 'resolveFamdbPath\(\)' "$script should define resolveFamdbPath()"
    assert_contains "$script" 'find "\$rm_prefix/share" -maxdepth 1 -type d -name "famdb-\*"' "$script should detect standalone famdb-* installs"
    assert_contains "$script" 'famdb\.py -i "\$libpath"' "$script should call famdb.py with resolved libpath"
    assert_contains "$script" 'touch "\$\{library_path\}/\.earlgrey\.config\.complete"' "$script should touch a concrete completion marker path"
done

assert_contains "scripts/LTR_FINDER_parallel" '^#!/usr/bin/env perl' 'LTR_FINDER_parallel should use env perl'
assert_contains "conda/meta.yaml" '^\s*-\s+blast\s*=2\.17\.0$' 'conda recipe should pin blast to 2.17.0'

echo "All regression checks passed."
