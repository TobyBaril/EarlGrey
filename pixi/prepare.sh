#!/usr/bin/env bash
#
# prepare.sh
#
# Applies the remaining source fixups from conda/build.sh that are needed to run
# EarlGrey straight from this repo under pixi (the env provides only the
# dependencies, not the packaged earlgrey).  Run alongside ./fix-scriptdir.sh
# (see the `setup` task in pixi.toml).  Idempotent.
#
# Steps (cf. conda/build.sh):
#   * Drop RepeatClassifier's -pa flag in scripts/TEstrainer/TEstrainer.
#     Newer RepeatModeler's RepeatClassifier no longer accepts -pa.
#   * Extract the bundled tRNA database that LTR_FINDER needs, and make the
#     ltr_finder binary executable.
#
# Skipped vs conda/build.sh (not applicable here):
#   * CONDA_DEFAULT_ENV name check removal - already absent from this source.
#   * `sed -i` -> `sed -i.bak` and sa-ssr `-t ${THREADS}` removal - macOS-only;
#     this env is linux-64.
#   * SCRIPT_DIR/STRAIN_SCRIPTS rewrite - handled by ./fix-scriptdir.sh.
#   * PERL5LIB activate/deactivate hooks - handled by [activation.env] in pixi.toml.
#   * SA-SSR build - handled by the `build-sassr` task.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

# --- RepeatClassifier -pa removal -------------------------------------------
TESTRAINER="$REPO/scripts/TEstrainer/TEstrainer"
if [ -f "$TESTRAINER" ]; then
    if grep -qF 'RepeatClassifier -pa ${THREADS} ' "$TESTRAINER"; then
        sed -i 's/RepeatClassifier -pa ${THREADS} /RepeatClassifier /' "$TESTRAINER"
        echo "patched: scripts/TEstrainer/TEstrainer (removed RepeatClassifier -pa)"
    else
        echo "already patched: scripts/TEstrainer/TEstrainer"
    fi
else
    echo "skip (missing): scripts/TEstrainer/TEstrainer" >&2
fi

# --- tRNAdb extraction for LTR_FINDER ---------------------------------------
LTRDIR="$REPO/scripts/bin/LTR_FINDER.x86_64-1.0.7"
if [ -f "$LTRDIR/tRNAdb.tar.gz" ]; then
    if [ -d "$LTRDIR/tRNAdb" ]; then
        echo "already extracted: tRNAdb"
    else
        tar -xzf "$LTRDIR/tRNAdb.tar.gz" -C "$LTRDIR"
        echo "extracted: scripts/bin/LTR_FINDER.x86_64-1.0.7/tRNAdb"
    fi
else
    echo "skip (missing): $LTRDIR/tRNAdb.tar.gz" >&2
fi
[ -f "$LTRDIR/ltr_finder" ] && chmod +x "$LTRDIR/ltr_finder" && echo "chmod +x ltr_finder"
