# Changes

Summary of additions and changes on this branch (`hyphaltip_explore_claudeImprove`)
compared to `main`.

- Added a pixi environment in pixi folder to support direct build and custom pkg load

## Bug fixes — RepeatCraft LTR merging (`scripts/repeatCraft/helper/truemergeltrm.py`)

The LTR true-merge helper had several latent bugs that could silently corrupt or
drop output. These are now fixed:

- **Final merged cluster was dropped.** The post-loop flush of the last
  accumulated cluster was commented out, so any GFF that ended while still inside
  an LTR group lost its final merged feature. The flush is restored (guarded by
  `len(d["col"]) > 0`, mirroring `truemergetem.py`) and now applies the merged
  start/end/strand before printing.
- **Strand-selection logic was dead.** The "use the strand of the largest
  fragment" check compared against `d["end"]` *after* it had already been
  overwritten, reducing the condition to one that is never true on sorted input.
  Replaced with an explicit running-maximum (`d["size"]`) so the merged feature
  correctly adopts the strand of its largest fragment.
- **Interactive `input()` corrupted output / crashed pipelines.** `sys.stdout`
  is redirected to the output GFF, so the "Skip this line?" prompt was being
  written *into* the result file, and `input()` raised `EOFError` in
  non-interactive runs (SLURM/Docker). The prompt now goes to `stderr` and
  `EOFError` defaults to skipping the malformed line.
- **Attribute parsing hardened.** `attr2dict` now uses `split("=", 1)` so
  attribute values containing `=` no longer raise `ValueError`.

### Malformed-line handling (committed earlier on this branch)
- Blank lines and rows with fewer than 9 tab-separated columns are now detected
  and skipped cleanly instead of triggering an `IndexError` deep in parsing.

## Bug fixes — RepeatCraft ErrorManagement variant

- **`truemergeltrmErrorManagement.py` was unimportable.** It contained a
  `TabError` (mixed tabs/spaces and a malformed `try/if/else` block) that
  prevented the module from being imported at all. The block was repaired into a
  valid `len(col) < 9` guard, and the same three fixes above
  (final-cluster flush, strand selection, `attr2dict`) were applied so it stays
  consistent with the primary helper.
- **Import/caller mismatch fixed.** `scripts/repeatcraftErrorManagement.py` and
  `scripts/repeatCraft/repeatcraftErrorManagement.py` imported
  `truemergeltrmErrorManagement` but then called `truemergeltrm.trumergeLTR(...)`
  — a module that is never imported in those files, which would raise a
  `NameError` whenever an LTR library was present. The call now correctly uses
  `truemergeltrmErrorManagement.trumergeLTR(...)`.

## Bug fixes — driver script (`earlGrey`)

- **`RepeatMasker -pa` underflow guard.** `rmthreads=$(( ProcNum / 4 ))` evaluated
  to `0` whenever fewer than 4 threads were requested, producing an invalid
  `RepeatMasker -pa 0`. All four occurrences now clamp to a minimum of 1
  (`(( rmthreads < 1 )) && rmthreads=1`).

## Performance — divergence calculation (`scripts/divergenceCalc/divergence_calc.py`)

Runtime improvements to the per-copy Kimura-80 divergence step, which runs once
per TE annotation in the genome (often 10^5–10^6 rows). All changes were verified
to produce **byte-identical output** versus the previous implementation.

- **Batched sequence extraction (one `bedtools getfasta` per chunk).** Previously
  every GFF row spawned its own `bedtools getfasta`, which re-opens the genome and
  its `.fai` index on each call. Sequences are now extracted with a single
  strand-aware `getfasta` call per worker chunk; the GFF row index is carried in
  the BED name column and recovered from the FASTA header, so the mapping is robust
  to record reordering and to header-format differences between bedtools versions.
  A whole-batch failure falls back to marking each affected row `NA`, matching the
  original per-row error handling.
- **Load balancing across workers.** The input was split into exactly one chunk
  per core, so `imap_unordered` could not rebalance: a worker handed a chunk full
  of slow or timeout-prone alignments became the long pole while the rest idled.
  The data is now split into `cores × chunk_factor` chunks (new `-c/--chunk_factor`
  option, default 4), so a worker that finishes early picks up more work. Capped
  at one chunk per row.
- **`Kimura80` lookup table.** The base-pair classifier now does a single O(1)
  dict lookup per aligned base (`_BASE_PAIR_CATEGORY`, mapping each ordered pair to
  `"MATCH"` / `"TRANSITION"` / `"TRANSVERSION"`) instead of three separate
  list-membership scans per base. The table is built once at import.
- **`Kimura80` math-domain crash fixed.** For saturated (highly divergent)
  alignments the Kimura-2-parameter terms `(1 − 2p − q)` and `(1 − 2q)` go to or
  below zero, where `log()`/`sqrt()` raised a `math domain error` and killed the
  worker. These cases are now detected and reported as `NA` (consistent with the
  unalignable case), so a single bad pair can no longer abort a chunk. Verified on
  100,000 random pairs: all 504 previously-crashing inputs now return `NA`, and all
  99,496 in-domain results are bit-for-bit unchanged.
- **Removed a bare `except:`** in the extraction path (now `except Exception`),
  and dropped the now-redundant per-row pybedtools cleanup counter.

### Benchmark

Measured in the repo `pixi/` environment (bedtools 2.31.1, samtools 1.23.1,
EMBOSS matcher 6.6.0).

**Batched sequence extraction.**

| Test | Original (per-row) | New (batched) | Speedup |
|------|--------------------|---------------|---------|
| Sequence extraction only, 600 features (600 → 1 subprocess spawns) | 3.13 s | 0.058 s | **~54×** |
| Full divergence run, 600 rows (synthetic 8 KB genome) | 10.27 s | 9.48 s | ~1.08× |

The full-run figure understates the gain: on the tiny synthetic genome the
per-copy `matcher` alignments (one subprocess per row, unchanged here) dominate
wall-clock, and `getfasta` is cheap because the `.fai` is trivially small. The
isolated extraction benchmark reflects the actual optimization, and its impact
grows with genome size (larger `.fai` to re-read per call) and with the number of
annotations. Correctness was confirmed identical at both 4 and 600 rows, across
plus/minus strands and multiple scaffolds.

**Load balancing (`--chunk_factor`).** Modelled with the real scheduling
machinery (`np.array_split` → `Pool(maxtasksperchild=1)` → `imap_unordered`),
2000 rows on 8 cores with a contiguous block of slow alignments (the worst case
for one-chunk-per-core). Total work 11.36 s; ideal balance ≈ 1.42 s.

| Chunking | Chunks | Wall | Speedup vs old |
|----------|--------|------|----------------|
| Old: chunks = cores | 8 | 9.64 s | 1.0× |
| New: chunks = cores × 4 (default) | 32 | 2.54 s | **3.8×** |
| New: chunks = cores × 8 | 64 | 1.47 s | **6.6×** |

The benefit appears only when per-repeat alignment times are uneven (which is
common — satellite-rich loci and hard-to-align families cluster in coordinate
order and dominate the `matcher` timeout budget). When work is uniform the change
is neutral.

## Tooling / repository hygiene

- **Added `modules/trf411.linux64`** — TRF 4.11 built from
  https://github.com/hyphaltip/TRF-64bit, which fixes a buffer overflow on long
  path names and applies aglabx's 64-bit integers for long chromosomes.
- **Removed tracked compiled artifacts.** 40 committed
  `scripts/repeatCraft/helper/__pycache__/*.pyc` files were removed from version
  control.
- **Added `.gitignore`** covering `__pycache__/`, `*.py[cod]`, and tool caches
  (`.ruff_cache/`, `.mypy_cache/`, `.pytest_cache/`) so build artifacts no longer
  get tracked.
- **Added `.vscode/settings.json`** enforcing tab indentation (this is a
  tab-indented Python project): `insertSpaces: false`, `detectIndentation:
  false`, `tabSize: 4`, applied globally and for `[python]`.
