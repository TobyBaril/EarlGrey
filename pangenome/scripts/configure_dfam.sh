#!/usr/bin/env bash
set -euo pipefail

FAMDB_DIR="$1"
RM_ROOT="$2"

echo "Configuring Dfam 3.9 in: $FAMDB_DIR"
cd "$FAMDB_DIR"

echo "Downloading Dfam partitions..."
curl -o 'dfam39_full.#1.h5.gz' \
  'https://dfam.org/releases/current/families/FamDB/dfam39_full.[0-16].h5.gz'

echo "Decompressing..."
gunzip -f *.gz

echo "Backing up min_init.0.h5"
if [[ -f min_init.0.h5 ]]; then
    mv min_init.0.h5 min_init.0.h5.bak
fi

echo "Running RepeatMasker configure..."
cd "$RM_ROOT/share/RepeatMasker"

perl ./configure \
    -libdir "$RM_ROOT/share/RepeatMasker/Libraries" \
    -trf_prgm "$RM_ROOT/bin/trf" \
    -rmblast_dir "$RM_ROOT/bin" \
    -hmmer_dir "$RM_ROOT/bin" \
    -abblast_dir "$RM_ROOT/bin" \
    -crossmatch_dir "$RM_ROOT/bin" \
    -default_search_engine rmblast

echo "Marking configuration complete"
touch "$FAMDB_DIR/.earlgrey.config.complete"

echo "Dfam 3.9 configuration complete."
