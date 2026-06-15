#!/usr/bin/env bash
#
# install.sh
#
# Installs EarlGrey *into the active pixi/conda environment*, mirroring the
# layout that conda/build.sh produces:
#
#   * stage the whole package under $CONDA_PREFIX/share/earlgrey-<version>
#     (the equivalent of conda's
#      ${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}-${PKG_BUILDNUM}),
#   * rewrite the baked-in SCRIPT_DIR / STRAIN_SCRIPTS to that location,
#   * apply the remaining source fixups (RepeatClassifier -pa, tRNAdb), and
#   * symlink the three earlGrey* executables into $CONDA_PREFIX/bin so they
#     are on PATH inside the env.
#
# SA-SSR (the `build-sassr` task) and the PERL5LIB activate/deactivate hooks
# ([activation.env] in pixi.toml) complete the install.
#
# Idempotent: re-running re-stages a fresh copy.
#
# Usage:   pixi run install     (preferred; also builds SA-SSR)
#          ./install.sh         (staging + symlinks only)
#
set -euo pipefail

: "${CONDA_PREFIX:?must be run inside an activated pixi/conda env (CONDA_PREFIX is unset)}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

# Package version from pixi.toml [workspace].
VERSION="$(sed -n 's/^[[:space:]]*version[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' "$HERE/pixi.toml" | head -n1)"
VERSION="${VERSION:-0.0.0}"

PKG_NAME="earlgrey"
PACKAGE_HOME="$CONDA_PREFIX/share/${PKG_NAME}-${VERSION}"
SCRIPT_DIR="$PACKAGE_HOME/scripts"

echo "Installing EarlGrey $VERSION -> $PACKAGE_HOME"

# --- stage the package ------------------------------------------------------
# Fresh copy each run (drops stale files from a previous version).
mkdir -p "$CONDA_PREFIX/bin"
rm -rf "$PACKAGE_HOME"
mkdir -p "$PACKAGE_HOME"
# Copy the repo, skipping VCS metadata and the pixi env / build caches.
tar -C "$REPO" \
    --exclude='./.git' \
    --exclude='./pixi/.pixi' \
    --exclude='./.ruff_cache' \
    --exclude='./pixi/.ruff_cache' \
    --exclude='./.trunk' \
    -cf - . | tar -C "$PACKAGE_HOME" -xf -

# --- rewrite SCRIPT_DIR / STRAIN_SCRIPTS to the installed location -----------
# (cf. conda/build.sh; matches assignment lines only, since references use
#  ${SCRIPT_DIR}, never SCRIPT_DIR= .)
for exe in earlGrey earlGreyAnnotationOnly earlGreyLibConstruct; do
    [ -f "$PACKAGE_HOME/$exe" ] && \
        sed -i "s|SCRIPT_DIR=.*|SCRIPT_DIR=${SCRIPT_DIR}|g" "$PACKAGE_HOME/$exe"
done
for helper in "$SCRIPT_DIR"/rcMergeRepeat* "$SCRIPT_DIR/headSwap.sh" "$SCRIPT_DIR/autoPie.sh"; do
    [ -f "$helper" ] && sed -i "s|SCRIPT_DIR=.*|SCRIPT_DIR=${SCRIPT_DIR}|g" "$helper"
done
STRAINER_WRAP="$SCRIPT_DIR/TEstrainer/TEstrainer_for_earlGrey.sh"
[ -f "$STRAINER_WRAP" ] && \
    sed -i "s|STRAIN_SCRIPTS=.*|STRAIN_SCRIPTS=${SCRIPT_DIR}/TEstrainer/scripts/|g" "$STRAINER_WRAP"

# --- remaining source fixups (cf. conda/build.sh, Linux-only subset) ---------
# RepeatClassifier -pa removal (newer RepeatModeler no longer accepts -pa).
TESTRAINER="$SCRIPT_DIR/TEstrainer/TEstrainer"
[ -f "$TESTRAINER" ] && \
    sed -i 's/RepeatClassifier -pa ${THREADS} /RepeatClassifier /' "$TESTRAINER"

# Extract the tRNA database LTR_FINDER needs.
LTRDIR="$SCRIPT_DIR/bin/LTR_FINDER.x86_64-1.0.7"
if [ -f "$LTRDIR/tRNAdb.tar.gz" ] && [ ! -d "$LTRDIR/tRNAdb" ]; then
    tar -xzf "$LTRDIR/tRNAdb.tar.gz" -C "$LTRDIR"
fi
rm -f "$LTRDIR/tRNAdb.tar.gz"

# --- permissions (cf. conda/build.sh) ---------------------------------------
chmod +x "$PACKAGE_HOME"/earlGrey "$PACKAGE_HOME"/earlGreyAnnotationOnly \
         "$PACKAGE_HOME"/earlGreyLibConstruct "$STRAINER_WRAP" 2>/dev/null || true
chmod +x "$SCRIPT_DIR"/* > /dev/null 2>&1 || true
[ -f "$LTRDIR/ltr_finder" ] && chmod +x "$LTRDIR/ltr_finder"
[ -d "$SCRIPT_DIR/repeatCraft/example" ] && chmod a+w "$SCRIPT_DIR/repeatCraft/example"

# --- precompile Python bytecode caches --------------------------------------
# Generate __pycache__/*.pyc now, as the install user, so the staged package
# under $PACKAGE_HOME (often read-only for those who later *run* EarlGrey) does
# not need to write bytecode at runtime. This benefits the helper modules that
# repeatcraft.py imports (fuseltr, fusetem, filtershortm, ...); entry scripts
# run as `python foo.py` never use a cache for themselves. `|| true`: some
# vendored sources may not compile cleanly and must not abort the install.
python -m compileall -q -j 0 "$PACKAGE_HOME" 2>/dev/null || true

# --- expose executables on PATH ---------------------------------------------
for exe in earlGrey earlGreyAnnotationOnly earlGreyLibConstruct; do
    ln -sf "$PACKAGE_HOME/$exe" "$CONDA_PREFIX/bin/$exe"
    echo "linked: \$CONDA_PREFIX/bin/$exe -> $PACKAGE_HOME/$exe"
done

echo "Done. 'earlGrey' is now on PATH inside the env."
