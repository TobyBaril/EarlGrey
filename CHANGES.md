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

- **Whitespace in `-s`/`-o` rejected up front.** A species name or output
  directory containing whitespace word-split the many unquoted
  `cd ${OUTDIR}/${species}_...` calls into multiple arguments, producing the
  cryptic `cd: too many arguments` failure mid-run. `Checks()` now validates both
  `$species` and `$directory` against `*[[:space:]]*` and exits with a clear
  message before any directory is created. Applied to all three entry points
  (`earlGrey`, `earlGreyAnnotationOnly`, `earlGreyLibConstruct`).

- **`looseMerge` directory creation is now idempotent.** `mergeRep()` ("Defragmenting
  Repeats") created its working directory with a bare
  `mkdir ${OUTDIR}/${species}_mergedRepeats/looseMerge`, which aborts the stage with
  `mkdir: cannot create directory ... File exists` on any rerun/resume where the
  directory already exists. Changed to `mkdir -p` (and quoted the path) so the stage
  re-enters cleanly. Applied to `earlGrey` and `earlGreyAnnotationOnly`; these were
  the only two `mkdir` calls in the codebase still missing `-p`.

- **Empty "latest run directory" globs no longer build garbage paths.** The
  `latestFile="$(realpath $(ls -td -- .../*/ | head -n 1) || true)/..."` idiom
  silently produced a bogus root-level path (e.g. `/<species>-families.fa.strained`)
  when the glob matched nothing — `|| true` swallowed the `realpath` error and the
  run continued with a broken path. Each site (TEstrainer in `strainer`/
  `strainerResume`, HELIANO resume) now captures the directory first, errors out
  with `exit 2` if none is found, and quotes the `realpath` argument; `ls` stderr
  is redirected to `/dev/null` so its "No such file" noise can no longer leak into
  the path. Applied across `earlGrey`, `earlGreyAnnotationOnly`, and
  `earlGreyLibConstruct`.

## New feature — checksum-validated resume (`earlGrey`)

- **`-v yes/no` (default `no`).** Re-running into an existing output directory has
  always skipped stages whose outputs are present, but that skip was based on
  output *existence* only — reusing a directory with a **different** input genome
  could silently reuse stale results. The new `-v yes` enables checksum-validated
  resume: when a stage finishes, the sha256 of the input genome (`$genOrig`) is
  recorded under `<output>/<species>_EarlGrey/.earlGrey_stamps/<stage>.sha256`, and
  on restart each stage (`initialMask`, `deNovo`, `strainer`, `novoMask`,
  `heliano`, `mergeRep`) is rerun if its recorded hash is missing or no longer
  matches the current genome. The genome is hashed once per run (cached).
- **Off by default and no-op when off.** The three helpers (`genomeHash`,
  `stageStale`, `stageStamp`) short-circuit when `-v` is not `yes`, so each stage
  guard reduces to the original `[ ! -f OUT ]` test and no stamp files are written.
  With the flag off the pipeline behaves byte-for-byte as before, keeping the fork
  aligned with upstream behaviour.
- **Caveats.** Enabling `-v yes` on a directory built without it forces one full
  rerun (no stamps yet). A stage judged stale re-executes over its existing stage
  directory and does not auto-clear prior outputs, so for a changed genome a fresh
  output directory remains the safest choice.

## Bug fixes — TEstrainer directory selection (`earlGrey`, `earlGreyLibConstruct`)

The `strainer()` and `strainerResume()` subprocesses locate the TEstrainer output
directory (`TS_<species>-families.fa_<timestamp>`) differently and had a few
fragile spots; both are now hardened and made consistent across the two driver
scripts (`earlGreyAnnotationOnly` has neither subprocess).

- **`strainerResume()` accepted a missing/empty directory.** The guard
  `[ ! $(echo "$strainDataDir" | wc -l) -eq 1 ]` treated the empty case as a
  single match (`echo ""` is one line), so a missing `TS_` directory passed the
  check and an empty `-d` was handed to TEstrainer (in `earlGreyLibConstruct` the
  `-z` check was absent entirely). The guard now distinguishes **zero** directories
  from **more than one**, with an accurate message for each, and counts matches via
  `grep -c .` instead of `wc -l` on a possibly-empty string.
- **Path was not fully resolved.** `strainerResume()` built
  `latestStrainDir=${OUTDIR}/${species}_strainer/${strainDataDir}`, where
  `find` returns a `./TS_...` prefix, producing a stray `/./` in the path. It now
  resolves the directory with `realpath` up front (while the cwd is still the
  strainer dir) so the downstream `latestFile` is a clean absolute path. The dead
  `[ -z "$latestStrainDir" ]` check (never true once the literal `$OUTDIR` prefix
  is present) was removed.
- **`strainer()` aligned with `strainerResume()`.** Directory discovery switched
  from `ls -td -- .../*/ | head -n 1` (which matches *any* subdirectory) to a
  `find -name "TS_${species}-families.fa_*" -printf '%T@\t%p\n' | sort -nr | head`
  selection that restricts to TEstrainer output dirs while preserving the
  "newest first" behaviour, then resolves with `realpath`.
- **Unsafe glob-into-`test` in the step 4 resume check (`earlGrey`).** The branch
  `elif [ -d ${OUTDIR}/${species}_strainer/TS_${species}-families.fa_* ]` fed a
  glob and unquoted variables to `test -d`: multiple matching directories expanded
  to several words and triggered `[: too many arguments` (the test then returned
  false and TEstrainer silently restarted from scratch instead of resuming), and
  whitespace in `$OUTDIR`/`$species` word-split the path. Replaced with a quoted
  `find -maxdepth 1 -type d -name "TS_${species}-families.fa_*"` whose result is
  tested with `[ -n ... ]`, so 0/1/many matches are all handled; the multi-dir case
  is left to `strainerResume()`, which already errors out explicitly.

## Bug fixes / cleanup — TEtrim.py (`scripts/TEstrainer/scripts/TEtrim.py`)

- **`TypeError` crash in debug mode.** Two `sys.exit(...)` messages concatenated
  the integer `args.minimum_seq` into a string, raising
  `TypeError: can only concatenate str (not "int") to str` whenever `-D TRUE` took
  the "too few sequences" branch. Both now wrap the value in `str()`.
- **Fragile Biopython version test.** `float(Bio.__version__) >= 1.79` (used twice)
  misparses versions like `1.100` (→ `1.1`). Replaced with a single module-level
  `_bio_modern` flag computed from an int-tuple parse of the version, de-duplicating
  the two `ungap("-")` vs `replace("-", "")` branches.
- **Shell-call path hardening.** The three `os.system(...)` invocations
  (blastn|uniq, mafft, final blastn) concatenated unquoted paths, so a space in the
  output directory or sequence name broke the command. Every interpolated path (and
  `args.threads`) is now wrapped in `shlex.quote(...)`, consistent with the
  whitespace-robustness work in the `earlGrey` wrappers; pipe/redirect behaviour is
  otherwise unchanged.
- **Minor tidy.** Removed the unused `import string`, simplified
  `align[1:len(align)]` to `align[1:]` and `if df.empty is True:` to `if df.empty:`,
  and stripped trailing whitespace.

## Bug fixes — RepeatCraft helpers (stderr messages)

- **Missing newlines on `stderr.write`.** The `"Updated LTR.gff with LTRgroup
  attribute to:"` status message in `scripts/repeatCraft/helper/fuseltr.py` and
  `repeatcraftHelper.py` was written without a trailing `\n`, so it ran into
  subsequent output. Both now terminate the line.
- **Progress counters left the terminal mid-line.** The `\rProgress:.../...`
  in-place counters in `extraFuseTEm.py`, `fusetem.py`, and both loops in
  `repeatcraftHelper.py` never emitted a final newline, so the next message was
  appended to the progress line. Each loop now writes a closing `\n`.

## Performance — TEstrainer per-family scripts (`scripts/TEstrainer/scripts/`)

These scripts are invoked once per TE family per iteration under GNU `parallel`
(thousands of invocations per run), so both per-invocation overhead and any
quadratic inner loops matter. The following are pure speed-ups with unchanged
output:

- **`initial_mafft_setup.py` — TRF parse was O(n²).** Each qualifying TRF line did
  `trf_df = pd.concat([trf_df, pd.DataFrame({...})])`, copying the whole frame
  every iteration. Coordinates are now collected into a Python list and the
  DataFrame is built once after the loop. Each line is also split once into
  `fields` instead of 2–4 times. As a side benefit the `Start`/`End` columns now
  build as native `int64` (the old empty-frame concat produced `object` dtype and
  a pandas FutureWarning) which is what the downstream `pr.PyRanges` call wants.
- **`initial_mafft_setup.py` — dead imports removed.** `string`, `statistics`, and
  `numpy as np` were imported but unused. The deliberate deferred-import block
  (pandas/pyranges loaded *after* the `file_check` calls, so missing-file
  invocations exit before paying the ~1 s import cost) was kept as-is.
- **`TEtrim.py` — `single_trim`/`con_maker` vectorised with NumPy.**
  `single_trim` rebuilt a new `MultipleSeqAlignment` by concatenating one column
  at a time (`good = good + aln_in[:,x:x+1]`), which is **O(length²)** and was the
  dominant cost on long TE+flank alignments. `con_maker` likewise looped per
  column in Python, re-scanning each column multiple times. Both now build the
  alignment into a single `(rows × cols)` char array once and do the column math
  vectorised: `single_trim` keeps columns with >1 non-gap base (case and record
  metadata preserved); `con_maker` emits a base only where >2 non-gap characters
  exist, with A/T/C/G ties (including all-zero columns) → `'n'` and otherwise the
  most frequent base in A,T,C,G order. Output is byte-identical to the previous
  implementation — verified with a throwaway harness over 4007 cases (a 4000-case
  random sweep plus targeted edge cases: ties, weak ≤2-non-gap columns, all-gap,
  all-N, mixed case, 2-row), all matching, then removed.
- **`Dfam_extractor.py` — library parsed twice.** The Dfam and non-Dfam splits
  each re-read and re-parsed the whole input FASTA. Now a single pass computes
  `is_dfam` once per record and writes to the matching handle (the partition is the
  exact logical complement of the original). Also fixed an adjacent guaranteed
  crash where the missing-directory branch called `os.mkdir(args.out_dir)` for a
  non-existent attribute (should be `args.directory`).
- **`splitter.py` / `trf_parser.py` — redundant `split()`.** `splitter.py` split
  each record name twice (now once into `base_name`); `trf_parser.py` re-split the
  header line after already computing `splitline` (now reuses `splitline[0]`).

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
- **Linear result aggregation.** `tmp_out_parser` concatenated each per-worker
  result file onto a growing DataFrame inside the loop, reallocating it every
  iteration (O(n²) in the number of chunks/rows). It now reads all files into a
  list and concatenates once (O(n)), which also avoids concatenating onto an empty
  DataFrame.
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

**Result aggregation (`tmp_out_parser`).** Concat-in-a-loop vs single concat,
scaling the number of per-worker result files (2000 rows each):

| Result files | Total rows | Old (O(n²)) | New (O(n)) | Speedup |
|--------------|-----------|-------------|------------|---------|
| 64 | 128,000 | 0.170 s | 0.102 s | 1.7× |
| 256 | 512,000 | 1.594 s | 0.371 s | 4.3× |
| 512 | 1,024,000 | 5.820 s | 0.758 s | **7.7×** |

The old timings grow ~4× per doubling (quadratic); the new ones ~2× (linear), so
the gap widens with genome size and thread count (both increase the chunk count).

## Tooling / repository hygiene

- **Added `modules/trf411.linux64`** — TRF 4.11 built from
  https://github.com/hyphaltip/TRF-64bit, which fixes a buffer overflow on long
  path names and applies aglabx's 64-bit integers for long chromosomes.
- **Removed tracked compiled artifacts.** 40 committed
  `scripts/repeatCraft/helper/__pycache__/*.pyc` files were removed from version
  control.
- **Precompile Python bytecode at install/build time.** Both `pixi/install.sh`
  and `conda/build.sh` now run `python -m compileall` over the staged package, so
  `__pycache__/*.pyc` caches are generated by the install user (the share
  directory is often read-only for those who later run EarlGrey). This benefits the
  helper modules `repeatcraft.py` imports; it is best-effort (`|| true`) so a
  vendored source that does not compile cannot abort the install/build.
- **Added `.gitignore`** covering `__pycache__/`, `*.py[cod]`, and tool caches
  (`.ruff_cache/`, `.mypy_cache/`, `.pytest_cache/`) so build artifacts no longer
  get tracked.
- **Added `.vscode/settings.json`** enforcing tab indentation (this is a
  tab-indented Python project): `insertSpaces: false`, `detectIndentation:
  false`, `tabSize: 4`, applied globally and for `[python]`.
