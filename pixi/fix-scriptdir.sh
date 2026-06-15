#!/usr/bin/env bash
#
# fix-scriptdir.sh
#
# EarlGrey's executables and helper scripts ship with the upstream developer's
# absolute path baked in, e.g.
#
#     SCRIPT_DIR=/data/toby/EarlGrey/scripts/
#     STRAIN_SCRIPTS=/data/toby/EarlGrey/scripts/TEstrainer/scripts/
#
# The conda build (conda/build.sh) rewrites these to the install prefix at build
# time. When running EarlGrey straight from this repo via pixi (the env provides
# only the dependencies, not the earlgrey package), those lines still point at
# the developer's machine and the pipeline fails on the first helper call.
#
# This script rewrites every such line to locate the scripts directory relative
# to the running script itself (resolving symlinks with readlink -f), so the
# code is correct no matter where the repo lives - no path is hardcoded. It is
# idempotent: re-running it (e.g. after pulling from upstream) is safe.
#
# Usage:   pixi/fix-scriptdir.sh
#          (run from anywhere; it locates the repo root from its own path)
#
set -euo pipefail

# Repo root = parent of the dir holding this script (pixi/).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

# Self-locating directory of the *running* script, symlink-resolved.
# Pre-escaped for use on the right-hand side of sed (the && contains &).
SELF_DIR='$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" \&\& pwd)'

# Marker used to detect an already-patched file.
MARKER='readlink -f "${BASH_SOURCE[0]}"'

# patch_var <file-relative-to-repo> <VAR> <suffix-after-SELF_DIR>
patch_var() {
    local rel="$1" var="$2" suffix="$3"
    local f="$REPO/$rel"

    if [ ! -f "$f" ]; then
        echo "skip (missing): $rel" >&2
        return
    fi
    if ! grep -qE "^[[:space:]]*${var}=" "$f"; then
        echo "WARNING: no ${var}= line in $rel (upstream layout changed?)" >&2
        return
    fi
    if grep -qF "$MARKER" "$f"; then
        echo "already patched: $rel"
        return
    fi
    sed -i -E "s|^([[:space:]]*)${var}=.*|\1${var}=\"${SELF_DIR}${suffix}\"|" "$f"
    echo "patched: $rel"
}

# Top-level executables live at the repo root; helpers are in ./scripts
for exe in earlGrey earlGreyAnnotationOnly earlGreyLibConstruct; do
    patch_var "$exe" SCRIPT_DIR "/scripts"
done

# Helper scripts live inside ./scripts, so SCRIPT_DIR is their own directory
for helper in \
    scripts/rcMergeRepeats \
    scripts/rcMergeRepeatsLoose \
    scripts/headSwap.sh \
    scripts/autoPie.sh; do
    patch_var "$helper" SCRIPT_DIR ""
done

# TEstrainer wrapper sits in ./scripts/TEstrainer, with its scripts one level down
patch_var scripts/TEstrainer/TEstrainer_for_earlGrey.sh STRAIN_SCRIPTS "/scripts"
