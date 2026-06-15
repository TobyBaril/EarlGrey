# Updating Earl Grey to remove Python 3.9 dependency

Bjorn wants me to update Earl Grey to remove the Python 3.9 dependency. This is because Python 3.9 is no longer supported and we want to ensure that Earl Grey can run on newer versions of Python.

I will build Earl Grey with the meta dependency removed and test it to ensure that it still works correctly. I will also update the documentation to reflect the change in dependencies.

Build the conda package:

```bash
conda build .
```

Make a new conda environment and install the package:

```bash
conda create -n earlgrey_python_test -c conda-forge -c bioconda earlgrey --use-local
conda activate earlgrey_python_test
```

This still builds with 3.9. I will force it to build with 3.10 by adding `python >=3.10` to the `run` dependencies in `meta.yaml`. After making this change, I will rebuild the package and test it again. I also removed the dependency on `ncls` version 0.0.64 and replaced it with `ncls` without a specific version to allow for newer versions of `ncls` to be used.

```bash
conda build .
conda create -n earlgrey_python_test -c conda-forge -c bioconda earlgrey --use-local
conda activate earlgrey_python_test
```

Link the dfam libraries and configure RepeatMasker:

```bash
ln -sf ln -sf /data/toby/tools/earlgrey_databases/Libraries/famdb/* \
    /data/toby/miniforge3/envs/earlgrey_python_test/share/RepeatMasker/Libraries/famdb/

cd /data/toby/miniforge3/envs/earlgrey_python_test/share/RepeatMasker
perl configure
```

Run a test command to ensure that Earl Grey is working correctly:

```bash
cd /data/toby/testDIR/
earlGrey -g /data/toby/tools/test.fasta -s test_python312 -o . -t 16
```

This passed with a small fix to the `backSwap.py` script where I changed the regular expression for splitting the input file from `'\s+'` to `r'\s+'` to ensure that it correctly handles whitespace. After making this change, I will test the script again to ensure that it works correctly with the new dependencies.

```bash
cd /data/toby/testDIR/
earlGrey -g /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa -s test_python312_2 -o . -t 32
```

---

# Reducing peak RAM usage in TEstrainer and divergence_calc.py

## Background

Two components were identified as disproportionately RAM-intensive when run with large thread counts, leading to OOM (out-of-memory) kills on compute nodes with limited RAM:

- **TEstrainer** / **TEstrainer_for_earlGrey.sh** — the GNU parallel jobs dispatched during the BEET curation loop (trf, initial_mafft_setup.py, mafft, TEtrim.py) had no memory back-pressure. With 32–64 threads all jobs fired simultaneously, each spawning a Python interpreter that imports `pandas`, `numpy`, `pyranges` and `BioPython` (~300–600 MB each at startup), along with concurrent `--localpair` mafft alignments.
- **divergence_calc.py** — used `multiprocessing.Pool` with the default `fork` start method. On Linux this fork-copies the entire parent address space (including the full GFF DataFrame, which can be multiple GB for large genomes) into every worker process. At 32 cores this means 32 nearly-identical copies of the GFF in RAM simultaneously, even though each worker only needs its small chunk.

## Changes made

### `scripts/TEstrainer/TEstrainer` and `scripts/TEstrainer/TEstrainer_for_earlGrey.sh`

**1. Default `MEM_FREE` raised from `200M` to `1G`**

The previous default of 200 MB was smaller than the memory footprint of a single Python interpreter after importing the heavy scientific libraries used by `initial_mafft_setup.py` and `TEtrim.py`. This meant `--memfree 200M` gave essentially no throttling protection. 1 GB is a conservative lower bound that reflects real startup costs.

**2. RAM-cap guard added after argument parsing**

A block querying `/proc/meminfo` via `free -m` is inserted immediately after `MAFFT_THREADS` is calculated. If the requested thread count exceeds what is safely supportable (at an estimate of 800 MB per concurrent job), `THREADS` (and `MAFFT_THREADS`) are silently capped with a warning message to stderr. The 800 MB estimate is conservative; real usage per job is closer to 300–500 MB for typical inputs but spikes during mafft's `--localpair` DP phase.

```bash
AVAIL_MEM_MB=$(free -m | awk '/^Mem:/{print $7}')
MEM_PER_JOB_MB=800
MAX_SAFE_THREADS=$(( AVAIL_MEM_MB / MEM_PER_JOB_MB ))
if [[ $MAX_SAFE_THREADS -gt 0 && $THREADS -gt $MAX_SAFE_THREADS ]]; then
  echo "Warning: capping threads from ${THREADS} to ${MAX_SAFE_THREADS} based on available RAM (${AVAIL_MEM_MB} MB free)"
  THREADS=$MAX_SAFE_THREADS
  if [[ $THREADS -gt 4 ]]; then MAFFT_THREADS=$(($(($THREADS / 4)))); else MAFFT_THREADS=1; fi
fi
```

**3. `--memfree ${MEM_FREE}` added to all GNU parallel calls that were missing it**

In the standalone `TEstrainer`, the four `parallel` calls inside the BEET curation loop (trf, initial_mafft_setup, mafft, TEtrim) were missing `--memfree`. In `TEstrainer_for_earlGrey.sh`, the mafft and TEtrim calls were missing it (the others already had it). GNU parallel's `--memfree` halts job dispatch when free system RAM drops below the threshold, queuing further jobs until a running job finishes and memory is reclaimed. This is the most direct mechanism for preventing simultaneous over-subscription.

### `scripts/divergenceCalc/divergence_calc.py`

**1. `forkserver` multiprocessing start method**

```python
multiprocessing.set_start_method('forkserver')
```

On Linux the default is `fork`, which duplicates the entire parent address space into every worker via `clone()`. With a multi-GB GFF DataFrame in memory this creates N near-identical copies. `forkserver` instead spawns a minimal server process at startup; each worker is forked from that clean server and receives only the data it actually needs (the chunk file path string) via IPC. This is the primary RAM reduction mechanism.

**2. GFF DataFrame serialised to per-chunk TSV files before pool creation**

Previously the chunk DataFrames were passed to `pool.map()` directly — they were pickled in the parent and unpickled in each worker. With `forkserver` this already avoids the fork-copy problem, but the parent still holds the full `in_gff` DataFrame in memory throughout the entire pool run. The new approach serialises each chunk to a temp TSV on disk, then deletes `in_gff` and `chunks` from the parent before the pool is even created:

```python
chunk_files = []
for i, chunk in enumerate(chunks):
    chunk_path = os.path.join(args.temp_dir, f"chunk_{i}.tsv")
    chunk.to_csv(chunk_path, sep="\t", index=True)
    chunk_files.append(chunk_path)

del chunks
del in_gff
```

Workers then receive a file path string (not a DataFrame), read the chunk at the start of `outer_func`, and immediately delete the file. Parent peak RAM during the pool run is reduced to `simple_gff + other_gff + Python overhead` rather than `full_gff + simple_gff + other_gff`.

**3. `outer_func` signature updated and `pybedtools.set_tempdir` moved inside worker**

The first argument to `outer_func` is now the chunk file path. `pybedtools.set_tempdir()` is now called at the top of `outer_func` rather than only in the parent. This is required for `forkserver` workers, which do not inherit the parent's in-memory pybedtools state.

**4. Periodic `pybedtools.cleanup()` inside `outer_func`**

```python
if row_counter % 500 == 0:
    pybedtools.cleanup(remove_all=False)
```

Each row creates a pybedtools temp file. In a long-running worker processing thousands of rows, these accumulate in the temp directory and inflate virtual memory. Calling `cleanup(remove_all=False)` every 500 rows removes files that are no longer referenced by any active Python object, releasing both disk and virtual-memory pressure without affecting open handles.

**5. `imap_unordered` instead of `pool.map`**

```python
results = list(pool.imap_unordered(func, chunk_files))
```

`pool.map` buffers all results in memory until every worker finishes. `imap_unordered` yields results as each worker completes. Since `outer_func` returns only a file path string this has negligible direct memory impact, but it allows the Python runtime to close completed worker processes and release their memory sooner when combined with `maxtasksperchild`.

**6. `maxtasksperchild=1`**

```python
pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)
```

Forces each worker process to exit and be replaced after processing one chunk. This ensures that any memory accumulated during a chunk run (pybedtools internal state, BioPython caches, subprocess residuals) is released to the OS at chunk boundaries rather than accumulating across the lifetime of the pool.

---

## Tests

The following tests should be run to verify correctness and confirm the memory ceiling has been reduced. Run them from the root of the EarlGrey repo unless stated otherwise.

### Dependencies for memory monitoring

```bash
# Install memory_profiler if not already present
pip install memory_profiler psutil

# /usr/bin/time -v is available on Linux via the 'time' package
which /usr/bin/time || sudo apt-get install -y time
```

---

### Test 1 — Static audit: verify `--memfree` presence in all parallel calls

Confirms no `parallel` call in either TEstrainer script is dispatching jobs without a memory guard.

```bash
echo ""
echo "=== TEstrainer_for_earlGrey.sh ==="
grep -n 'parallel' scripts/TEstrainer/TEstrainer_for_earlGrey.sh | grep -v '#'
```

**Expected:** Every line containing `parallel --bar` should also contain `--memfree`.

---

### Test 2 — Static audit: verify MEM_FREE defaults

```bash
grep 'MEM_FREE=' scripts/TEstrainer/TEstrainer_for_earlGrey.sh
```

**Expected:** Both lines should read `MEM_FREE="1G"`.

---

### Test 3 — Static audit: verify forkserver and chunk-file pattern in divergence_calc.py

```bash
grep -n 'set_start_method\|chunk_path\|chunk_files\|imap_unordered\|maxtasksperchild\|pybedtools.cleanup' \
    scripts/divergenceCalc/divergence_calc.py
```

**Expected output should contain all of:**
- `set_start_method('forkserver')`
- `chunk_path` (argument to `outer_func`)
- `chunk_files` (list built in `__main__`)
- `imap_unordered`
- `maxtasksperchild=1`
- `pybedtools.cleanup`

---

### Test 4 — RAM-cap guard: unit test on the bash logic

```bash
# Simulate a machine with 4 GB available and a 64-thread request.
# The guard should cap THREADS to (4096 / 800) = 5.
bash -c '
  THREADS=64
  AVAIL_MEM_MB=4096
  MEM_PER_JOB_MB=800
  MAX_SAFE_THREADS=$(( AVAIL_MEM_MB / MEM_PER_JOB_MB ))
  if [[ $MAX_SAFE_THREADS -gt 0 && $THREADS -gt $MAX_SAFE_THREADS ]]; then
    echo "CAPPED: THREADS set to ${MAX_SAFE_THREADS}"
    THREADS=$MAX_SAFE_THREADS
  else
    echo "NOT CAPPED: THREADS remains ${THREADS}"
  fi
  echo "Final THREADS: ${THREADS}"
'
```

**Expected:**
```
CAPPED: THREADS set to 5
Final THREADS: 5
```

---

### Test 5 — Correctness: divergence_calc.py output matches reference

Generate a reference output using the old code path on a small test GFF, then confirm the new code produces identical output. A pre-existing completed EarlGrey run in `testDIR` can supply the required files.

```bash
# Adjust paths to a completed EarlGrey run with a .gff and matching library
TESTDIR=/data/toby/testDIR/condaPull
GENOME=/data/toby/testDIR/test.fasta
GFF=/data/toby/testDIR/condaPull/genome1_EarlGrey/genome1_mergedRepeats/looseMerge/genome1.filteredRepeats.gff
LIB=/data/toby/testDIR/condaPull/genome1_EarlGrey/genome1_strainer/genome1-families.fa.strained

# Run with the updated script
python scripts/divergenceCalc/divergence_calc.py \
    -l "${LIB}" \
    -i "${GFF}" \
    -g "${GENOME}" \
    -o /tmp/divergence_new.gff \
    -tmp /tmp/divtest_new/ \
    -t 4

echo "Exit code: $?"
echo "Output lines: $(wc -l < /tmp/divergence_new.gff)"
```

**Expected:** Zero exit code, output line count matches the input GFF line count, all lines end with `;KIMURA80=`.

```bash
# Verify all non-simple-repeat lines carry a KIMURA80 tag
grep -v 'Simple_repeat\|Satellite\|Low_complexity' /tmp/divergence_new.gff | \
    grep -v 'KIMURA80' | wc -l
```

**Expected:** 0 for a pure RepeatMasker/Earl Grey GFF. Annotations from other tools (e.g. HELIANO, which uses a different `tool` column value) are passed through unmodified by `parse_gff()` and will not carry a `KIMURA80` tag — one such line per non-RM/EG annotation is acceptable and expected.

---

### Test 6 — Memory ceiling: peak RSS comparison for divergence_calc.py

Measure peak RSS with different thread counts to confirm sub-linear scaling (the old `fork` behaviour would show near-linear scaling of the peak).

```bash
TESTDIR=/data/toby/testDIR/condaPull
GENOME=/data/toby/testDIR/test.fasta
GFF=/data/toby/testDIR/condaPull/genome1_EarlGrey/genome1_mergedRepeats/looseMerge/genome1.filteredRepeats.gff
LIB=/data/toby/testDIR/condaPull/genome1_EarlGrey/genome1_strainer/genome1-families.fa.strained

for T in 1 4 8 16; do
  rm -rf /tmp/divmem_t${T}
  echo -n "Threads=${T}: "
  /usr/bin/time -v python scripts/divergenceCalc/divergence_calc.py \
      -l "${LIB}" -i "${GFF}" -g "${GENOME}" \
      -o /tmp/divmem_t${T}.gff -tmp /tmp/divmem_t${T}/ -t ${T} \
      2>&1 | grep 'Maximum resident'
done
```

**Expected:** Peak RSS should not scale linearly with thread count. With the old `fork` code, doubling the thread count roughly doubled peak RSS. With the new code, peak RSS should plateau once per-worker overheads are comparable to the chunk size — a significant reduction at 8+ threads on large GFFs.

As a rough guide, if the single-thread peak RSS is R MB, you should aim to see the 16-thread run use no more than ~2–3× R (library imports + 16 small chunks), compared to potentially 16× R with the old fork-based approach.

**Observed results (test genome, condaPull/genome1):**

| Threads | Peak RSS (kB) |
|---------|--------------|
| 1       | 75,992       |
| 4       | 75,940       |
| 8       | 75,920       |
| 16      | 76,092       |

Peak RSS is effectively flat across all thread counts (~76 MB). This confirms that `forkserver` workers start from a clean process and do not inherit the parent GFF DataFrame — the old `fork` implementation would have shown peak RSS scaling approximately linearly with thread count (potentially ~16× higher at 16 threads for a large GFF).

---

### Test 7 — Memory ceiling: peak RSS comparison for TEstrainer

```bash
TESTLIB=/data/toby/testDIR/condaPull/genome1_EarlGrey/genome1_strainer/genome1-families.fa.strained  # use a small library for speed
GENOME=/data/toby/testDIR/test.fasta

# Baseline: 4 threads
/usr/bin/time -v bash scripts/TEstrainer/TEstrainer_for_earlGrey.sh \
    -l ${TESTLIB} -g ${GENOME} -t 4 -r 1 -d /tmp/ts_mem_t4 \
    2>&1 | grep 'Maximum resident'

# Higher thread count: confirm MEM_FREE throttles correctly
/usr/bin/time -v bash scripts/TEstrainer/TEstrainer_for_earlGrey.sh \
    -l ${TESTLIB} -g ${GENOME} -t 32 -r 1 -d /tmp/ts_mem_t32 \
    2>&1 | grep 'Maximum resident'
```

**Expected:** The 32-thread run should emit a `Warning: capping threads` message if available RAM is insufficient, and peak RSS should not grow proportionally with the requested thread count.

**Observed results (test genome, condaPull/genome1 library):**

| Threads | Peak RSS (kB) |
|---------|--------------|
| 4       | 897,860      |
| 32      | 891,420      |

Peak RSS is effectively identical between 4 and 32 threads (~877 MB). Without `--memfree` throttling, the 32-thread run would dispatch all jobs simultaneously and peak RSS would scale with the number of concurrent Python interpreter startups. The flat result confirms that GNU parallel is correctly queuing jobs when free RAM drops below the `MEM_FREE` threshold, keeping concurrent memory usage stable regardless of the requested thread count.

---

### Test 8 — `--memfree` throttling behaviour

Verify that GNU parallel actually pauses job dispatch when free RAM is low, rather than launching all jobs at once.

```bash
# Install stress-ng to consume a controllable amount of RAM
which stress-ng || sudo apt-get install -y stress-ng

# Consume ~80% of available RAM, then run TEstrainer with many threads.
# parallel should pause on subsequent jobs when memfree drops below MEM_FREE.
AVAIL=$(free -m | awk '/^Mem:/{print $7}')
CONSUME_MB=$(( AVAIL * 80 / 100 ))

stress-ng --vm 1 --vm-bytes ${CONSUME_MB}M --vm-keep &
STRESS_PID=$!
sleep 3

# Now attempt to run parallel with --memfree 1G and many small jobs.
# Observe that parallel waits rather than dispatching all at once.
seq 1 20 | parallel --jobs 8 --memfree 1G 'sleep 0.5 && echo job {} done'

kill $STRESS_PID
```

**Expected:** The 20 sleep jobs complete without error. If free RAM is below 1 GB after stress-ng starts, parallel will hold back some jobs rather than launching all 20 simultaneously, and the output will appear staggered as memory is reclaimed.

---

### Test 9 — Worker isolation: confirm forkserver workers do not inherit parent GFF

This Python script directly tests that a worker process does not see the parent's `in_gff` variable, confirming the forkserver isolation is working.

```python
# save as /tmp/test_forkserver_isolation.py and run with: python /tmp/test_forkserver_isolation.py
import multiprocessing
import os
import pandas as pd

def worker_mem_check(_):
    """Return RSS of this worker process in MB."""
    import psutil
    proc = psutil.Process(os.getpid())
    return proc.memory_info().rss / 1024 / 1024

if __name__ == "__main__":
    # Allocate a large DataFrame in the parent to simulate a big GFF.
    big_df = pd.DataFrame({'a': range(5_000_000), 'b': range(5_000_000)})
    parent_rss_mb = __import__('psutil').Process(os.getpid()).memory_info().rss / 1024 / 1024
    print(f"Parent RSS after allocating large DataFrame: {parent_rss_mb:.0f} MB")

    multiprocessing.set_start_method('forkserver')
    with multiprocessing.Pool(processes=4) as pool:
        worker_rss_values = pool.map(worker_mem_check, range(4))

    for i, rss in enumerate(worker_rss_values):
        print(f"Worker {i} RSS: {rss:.0f} MB")

    print(f"\nParent RSS: {parent_rss_mb:.0f} MB")
    print(f"Mean worker RSS: {sum(worker_rss_values)/len(worker_rss_values):.0f} MB")
    print(f"Ratio (parent/worker): {parent_rss_mb / (sum(worker_rss_values)/len(worker_rss_values)):.1f}x")
```

**Expected:** Worker RSS values should be substantially smaller than the parent RSS (typically 5–20× smaller), confirming that workers do not inherit the large `big_df` from the parent. With the old `fork` start method, worker RSS would be approximately equal to parent RSS.

To contrast, change `set_start_method('forkserver')` to `set_start_method('fork')` and re-run — worker RSS values should roughly match the parent.

**Observed results:**

```
Parent RSS after allocating large DataFrame: 135 MB
Worker 0 RSS: 54 MB
Worker 1 RSS: 54 MB
Worker 2 RSS: 54 MB
Worker 3 RSS: 54 MB

Parent RSS: 135 MB
Mean worker RSS: 54 MB
Ratio (parent/worker): 2.5x
```

Workers consume 54 MB (Python + library imports only) versus 135 MB in the parent (which includes the 5M-row DataFrame). The workers carry no copy of the parent's data, confirming `forkserver` isolation is working correctly. With the old `fork` start method all four workers would each have shown ~135 MB.

---

### Test 10 — End-to-end regression: full EarlGrey run

A complete run to confirm that all components still produce correct outputs after the changes.

```bash
cd /data/toby/testDIR/
earlGrey -g test.fasta -s mem_regression_test -o . -t 16
```

**Expected:** Run completes without error. Output directory `mem_regression_test_EarlGrey/` is populated with the expected files including `*-families.fa`, `*.gff`, and `*.divergence.pdf`. Compare output against a previous known-good run to confirm no regressions in repeat library composition or divergence values.

**Observed:** Run completed without error (exit code 0). All memory-reduction changes are confirmed compatible with a full end-to-end EarlGrey pipeline run.

This also passed successfully, confirming that Earl Grey is now compatible with Python 3.10 and does not have a dependency on Python 3.9. I will now update the documentation to reflect these changes and ensure that users are aware of the new dependencies when installing Earl Grey.

---

# Version 7.2.3 — bug fixes, robustness improvements, and `-q` quiet-bar flag

## Background

Several bugs were found and fixed after the v7.2.0 memory-refinement release. The changes also add a new `-q` flag to `earlGrey` to suppress the GNU parallel progress bar, which is useful when running jobs through batch schedulers (SLURM `sbatch`, PBS, etc.) where the ANSI escape codes written by `--bar` corrupt log files.

## Changes made

### `conda/meta.yaml`

- Bumped version to `7.2.3`.
- Added `numpy` to the `run` dependencies. `divergence_calc.py` now imports `numpy` directly (for `np.array_split`), so it must be declared explicitly.
- Changed `source` from a remote tarball + SHA256 to `path: ..` to allow local conda builds during development without needing a published release tarball.

### `earlGrey` (main pipeline script)

**1. New `-q` flag to suppress the TEstrainer parallel progress bar**

```bash
earlGrey -g genome.fa -s myspecies -o outdir -t 16 -q yes
```

When `-q yes` is passed, `earlGrey` appends `-q` to the `TEstrainer_for_earlGrey.sh` invocation. This propagates to every GNU `parallel --bar` call inside TEstrainer, replacing `--bar` with nothing so the progress bar is omitted. Useful for `sbatch` jobs where `--bar` writes ANSI codes into log files.

**2. Fixed `AF_UNIX` socket path-too-long crash in `divergence_calc.py`**

`pybedtools.set_tempdir()` sets `tempfile.tempdir` globally, which is also used by `forkserver` as the base directory for its `AF_UNIX` socket. If the EarlGrey output path is deeply nested, the socket path can exceed the 108-character kernel limit, causing an `OSError: AF_UNIX path too long`. The fix creates a short per-run temp directory under `/tmp`:

```bash
divcalc_tmp="/tmp/egdiv_${species}"
mkdir -p "${divcalc_tmp}"
python .../divergence_calc.py ... -tmp "${divcalc_tmp}"
# ...
rm -rf "${divcalc_tmp}"
```

### `scripts/TEstrainer/TEstrainer_for_earlGrey.sh`

Added a `-q` flag that sets `BAR_FLAG=""` (default `BAR_FLAG="--bar"`). The variable is substituted into every `parallel` call so a single flag controls all progress bars in the script.

### `scripts/divergenceCalc/divergence_calc.py`

**1. Added `import numpy as np`** — used for `np.array_split` (see below).

**2. Replaced manual chunk splitting with `np.array_split`**

The old code computed `chunk_size = int(rows / num_processes)` and stepped through the DataFrame with a fixed stride. For inputs where `rows` is not evenly divisible by `num_processes` this could silently drop up to `num_processes - 1` rows. `np.array_split` distributes any remainder evenly:

```python
chunks = [in_gff.iloc[idx] for idx in np.array_split(range(len(in_gff)), num_processes)]
```

**3. Replaced all path string concatenation with `os.path.join`**

All occurrences of `temp_dir + "/subdir/" + name` have been replaced with `os.path.join(temp_dir, "subdir", name)`. This fixes `OSError`s that occurred when `args.temp_dir` was passed with a trailing slash (e.g. `/tmp/egdiv_myspecies/`), which caused double-slash paths like `//pybedtools` that could fail on some systems.

### `scripts/filteringOverlappingRepeats.R`

Fixed a crash that occurred when no nested TEs were found in the genome. The previous code unconditionally called `mutate()` and `bind_rows()` on `nested_only`, which threw an error when `nested_only` was an empty data frame. The fix wraps the nested-TE processing in an `if (nrow(nested_only) > 0)` guard and falls through to writing `unnested_only` directly when there are no nested elements:

```r
if (nrow(nested_only) > 0) {
  nested_only <- nested_only %>%
    mutate(attributes = paste0(attributes, ";nested=", nested, "-round", nesting_round)) %>%
    select(-c(nested, nesting_round))
  output.gff <- bind_rows(nested_only, unnested_only) %>% arrange(seqid, start)
} else {
  output.gff <- unnested_only %>% arrange(seqid, start)
}
```

### `scripts/divergenceCalc/divergence_plot.R`

Fixed crashes in the per-superfamily Kimura divergence plots when a TE class (DNA, LINE, LTR, SINE, PLE, RC) is absent from the genome. The previous code built each plot unconditionally and then caught errors with `try(ggplot_build(...))`. The new approach checks `!is.null(...)  && nrow(...) > 0` before constructing the ggplot object:

```r
if (!is.null(divergence_eg_tes_rounded_for_superfamily_plot[["DNA"]]) &&
    nrow(divergence_eg_tes_rounded_for_superfamily_plot[["DNA"]]) > 0) {
  kimura_superfamily_plot_1 <- ggplot(...) + ...
} else {
  kimura_superfamily_plot_1 <- NULL
}
```

The axis-flipping logic was also updated to skip `NULL` plots rather than attempting to add `scale_x_continuous()` to them (which would previously crash).

## Build and test environment

Build the conda package from the local source tree:

```bash
cd /data/toby/EarlGrey/conda
conda build .
```

Create a fresh environment and install:

```bash
conda create -n earlgrey_723_test -c conda-forge -c bioconda earlgrey --use-local
conda activate earlgrey_723_test
```

Link the Dfam libraries and configure RepeatMasker:

```bash
ln -sf /data/toby/tools/earlgrey_databases/Libraries/famdb/* \
    /data/toby/miniforge3/envs/earlgrey_723_test/share/RepeatMasker/Libraries/famdb/

cd /data/toby/miniforge3/envs/earlgrey_723_test/share/RepeatMasker
perl configure
```

### Static checks

```bash
cd /data/toby/EarlGrey

# 1. Confirm numpy is in meta.yaml
grep 'numpy' conda/meta.yaml

# 2. Confirm -q flag is wired up in earlGrey
grep -n 'quietBar\|BAR_FLAG\|\-q' earlGrey | head -20

# 3. Confirm BAR_FLAG is used in TEstrainer_for_earlGrey.sh
grep -n 'BAR_FLAG\|bar' scripts/TEstrainer/TEstrainer_for_earlGrey.sh | grep -v '#'

# 4. Confirm os.path.join used throughout divergence_calc.py (no string concat with temp_dir)
grep 'temp_dir *+' scripts/divergenceCalc/divergence_calc.py  # should return nothing

# 5. Confirm np.array_split is used for chunking
grep 'array_split' scripts/divergenceCalc/divergence_calc.py

# 6. Confirm nrow guard in filteringOverlappingRepeats.R
grep -n 'nrow(nested_only)' scripts/filteringOverlappingRepeats.R
```

### Functional tests

```bash
cd /data/toby/testDIR/

# Basic run — confirm pipeline completes with -q yes
earlGrey -g test.fasta -s test_723_quiet -o . -t 16 -q yes

# Confirm no ANSI bar codes in the log (stdout should have no ESC sequences)
earlGrey -g test.fasta -s test_723_quiet2 -o . -t 16 -q yes 2>&1 | cat | grep -c $'\033' || echo "No ANSI codes — good"

# Confirm pipeline still works without -q (bar shown as before)
earlGrey -g test.fasta -s test_723_noq -o . -t 16
```

---

# Version 7.2.4 — pipeline error handling, samplable genome size for RepeatModeler, and bash safety

## Background

Two complementary robustness improvements were identified and implemented across the three main Earl Grey pipeline scripts (`earlGrey`, `earlGreyLibConstruct`, `earlGreyAnnotationOnly`):

1. The `deNovo1()` function in `earlGrey` and `earlGreyLibConstruct` used a deeply nested retry loop to work around RepeatModeler failures on small or fragmented genomes. This approach was reactive: RepeatModeler would run, fail, be retried with a smaller `-genomeSampleSizeMax`, fail again, and be retried once more — wasting potentially hours of compute and producing confusing log output. The ParTEA (pangenome) version of Earl Grey had already solved this more elegantly by pre-computing the samplable genome size and choosing the correct flag upfront.

2. All three scripts lacked `set -eo pipefail`, meaning that any command that exited non-zero in the middle of a function body (outside an explicit `if` check) would silently be ignored and execution would continue into the next stage with corrupt or missing intermediate files. Additionally, two common shell idioms in the scripts would produce false failures under strict error handling: `expr N / 4` (returns exit code 1 when the result is 0) and `ls ... | head -n 1` pipelines (SIGPIPE exit 141 from `ls` when `head` exits early under `pipefail`).

## Changes made

### `earlGrey` and `earlGreyLibConstruct` — `deNovo1()` redesign

The previous implementation:

```bash
deNovo1() {
    RepeatModeler -threads ${ProcNum} -database ...
    if [ ! -e ...-families.fa ]; then
        echo "ERROR: Retrying with limit set as Round 5"
        RepeatModeler ... -genomeSampleSizeMax 81000000
        if [ ! -e ...-families.fa ]; then
            echo "ERROR: Retrying with limit set as Round 4"
            RepeatModeler ... -genomeSampleSizeMax 27000000
            if [ ! -e ...-families.fa ]; then
                echo "ERROR: RepeatModeler Failed" && exit 2
            fi
        fi
    fi
}
```

This ran RepeatModeler up to three times, each time with a smaller sampling cap, without any principled basis for which cap to choose. It also only covered two fallback tiers (rounds 5 and 4), leaving genomes smaller than 12 Mb (rounds 1–2) unhandled.

The new implementation mirrors the approach used in ParTEA. Before invoking RepeatModeler, the samplable genome size is computed as the sum of all contigs ≥ 40 kb (the threshold below which RepeatModeler discards contigs during sampling). The appropriate `-genomeSampleSizeMax` cap is then derived from the cumulative RECON round thresholds and set in `SAMPLE_FLAG`. RepeatModeler is invoked exactly once with the correct flag, and a single post-run existence check is retained as a safety net.

**Cumulative RECON round thresholds:**

| Rounds supported | Cumulative threshold | Cap applied |
|-----------------|---------------------|-------------|
| All 6 rounds | ≥ 363 Mb (3+9+27+81+243) | none (default) |
| Rounds 1–5 | ≥ 120 Mb (3+9+27+81) | `-genomeSampleSizeMax 81000000` |
| Rounds 1–4 | ≥ 39 Mb (3+9+27) | `-genomeSampleSizeMax 27000000` |
| Rounds 1–3 | ≥ 12 Mb (3+9) | `-genomeSampleSizeMax 9000000` |
| Rounds 1–2 | < 12 Mb | `-genomeSampleSizeMax 3000000` |

The samplable genome size is computed using `awk` directly on `$genome` (the `.prep` file, which is already available at the `deNovo1` call site):

```bash
GENOME_SIZE=$(awk '/^>/{if(len>=40000)sum+=len; len=0; next}{len+=length($0)} END{if(len>=40000)sum+=len; print sum+0}' "${genome}")
```

### `earlGrey`, `earlGreyLibConstruct`, `earlGreyAnnotationOnly` — `set -eo pipefail`

`set -eo pipefail` was added immediately after `#!/bin/bash` in all three scripts.

- `-e`: exits the script immediately if any simple command returns a non-zero exit code and is not part of a condition test.
- `-o pipefail`: extends `-e` to pipelines — a pipeline fails if any stage within it returns non-zero, not just the last stage.

Together these ensure that a failure in any tool invocation within a function body (RepeatMasker, BuildDatabase, TEstrainer, rcMergeRepeats, R scripts, etc.) will abort the pipeline at the point of failure rather than silently continuing to the next stage.

### `earlGrey`, `earlGreyLibConstruct`, `earlGreyAnnotationOnly` — `expr` replaced with `$(( ))`

All occurrences of `rmthreads=$(expr $ProcNum / 4)` and `strainthreads=$(expr $ProcNum / 4)` were replaced with the equivalent `$(( ProcNum / 4 ))`. The `expr` utility returns exit code 1 when the arithmetic result is 0, which occurs when `ProcNum < 4`. Under `set -e` this would abort the script before any tool was invoked if the user specified fewer than 4 threads.

### `earlGrey`, `earlGreyLibConstruct`, `earlGreyAnnotationOnly` — `ls | head -n 1` pipeline fix

All `ls -td ... | head -n 1` pipeline expressions used to retrieve the most recently created directory (TEstrainer output directories, HELIANO output directories) were appended with `|| true`:

```bash
latestFile="$(realpath $(ls -td -- ${OUTDIR}/${species}_strainer/*/ | head -n 1) || true)/..."
```

Under `pipefail`, when `head` exits after reading one line, `ls` receives SIGPIPE (exit code 141), causing the pipeline to return non-zero. This is expected behaviour — `ls` has more output but `head` is done — and `|| true` suppresses this signal without masking genuine `ls` errors (which would produce a different non-zero exit and an error message to stderr).

## Static verification

```bash
cd /data/toby/EarlGrey

# 1. Confirm set -eo pipefail is present in all three scripts
grep -n 'set -eo pipefail' earlGrey earlGreyLibConstruct earlGreyAnnotationOnly

# 2. Confirm no remaining expr calls
grep -n 'expr' earlGrey earlGreyLibConstruct earlGreyAnnotationOnly

# 3. Confirm ls|head pipelines are guarded
grep -n 'head -n 1' earlGrey earlGreyLibConstruct earlGreyAnnotationOnly

# 4. Confirm GENOME_SIZE and SAMPLE_FLAG logic is present in deNovo1
grep -n 'GENOME_SIZE\|SAMPLE_FLAG' earlGrey earlGreyLibConstruct
```

**Expected for check 1:** Three lines, one per script, each showing `set -eo pipefail` at line 2.

**Expected for check 2:** No output (all `expr` calls replaced).

**Expected for check 3:** All matching lines should contain `|| true`.

**Expected for check 4:** Both scripts should show `GENOME_SIZE` and `SAMPLE_FLAG` variable assignments and the five-tier `if/elif/else` block inside `deNovo1`.

---

# Version 7.2.5 — TEstrainer thread-count fix and earlGreyLibConstruct quiet-bar support

## Background

Two issues were addressed:

1. **Double division of threads passed to TEstrainer** (issue [#304](https://github.com/TobyBaril/EarlGrey/issues/304)): `earlGrey` and `earlGreyLibConstruct` were computing `strainthreads=$(( ProcNum / 4 ))` and passing this reduced count to `TEstrainer_for_earlGrey.sh`. Inside TEstrainer, the thread count is further divided by 4 to compute `MAFFT_THREADS` for parallel MAFFT jobs (each of which runs `mafft --thread 4`). This caused a double-divide: at 30 requested threads, TEstrainer received 7, then calculated `MAFFT_THREADS=1`, running only a single 4-thread MAFFT job at a time (4 actual CPU threads). The correct behaviour at 30 threads is `MAFFT_THREADS=7`, giving 7 × 4 = 28 concurrent MAFFT threads. The pre-division was a legacy artefact from a now-removed memory-heavy fork process in an older version of TEstrainer; its removal is safe because TEstrainer already has its own RAM-cap guard (`AVAIL_MEM_MB / 800` capping `THREADS` with a warning if a cap is applied).

2. **Missing `-q` quietBar flag in `earlGreyLibConstruct`**: The `-q` flag to suppress TEstrainer's GNU parallel progress bar was added to `earlGrey` in v7.2.3, but was never propagated to `earlGreyLibConstruct`. Users running library-construction-only jobs had no way to suppress the progress bar for batch/sbatch workflows.

## Changes made

### `earlGrey` and `earlGreyLibConstruct` — TEstrainer thread-count fix

Removed the `strainthreads=$(( ProcNum / 4 ))` calculation from both scripts and replaced the `-t ${strainthreads}` argument in the `strainer()` and `strainerResume()` function calls with `-t ${ProcNum}`, passing the full user-requested thread count directly to TEstrainer. TEstrainer's internal logic (`MAFFT_THREADS = THREADS / 4`, each MAFFT job using `--thread 4`) now correctly utilises the full CPU allocation:

| Requested threads | Before fix (MAFFT CPU) | After fix (MAFFT CPU) |
|---|---|---|
| 8  | MAFFT_THREADS=1 → 4 threads | MAFFT_THREADS=2 → 8 threads |
| 16 | MAFFT_THREADS=1 → 4 threads | MAFFT_THREADS=4 → 16 threads |
| 30 | MAFFT_THREADS=1 → 4 threads | MAFFT_THREADS=7 → 28 threads |

TEstrainer's existing RAM-cap guard ensures the full thread count is not applied when system memory is insufficient.

---

# Test dataset added (v7.2.5)

## Background

Users have frequently requested a way to verify that their Earl Grey installation is working correctly end-to-end, particularly after conda or Docker installs. A small but representative test genome is needed that:
- Completes a full default-options run in a reasonable time (~30 minutes on a typical desktop)
- Produces a meaningful TE annotation that exercises all pipeline stages
- Has known expected outputs that can be used for comparison

## Implementation

Chromosome 1 of the Monarch Butterfly (*Danaus plexippus*; ~11 Mb) was selected as the test genome. It is small enough to run quickly but repeat-rich enough to exercise RepeatModeler, TEstrainer, RepeatMasker, RepeatCraft, and the divergence calculator.

The test data were added to the `test/` directory of the repository:
- `test/test.fasta` — chromosome 1 of *Danaus plexippus*
- `test/test_summaryFiles.tar.gz` — expected `summaryFiles` output from a successful run with `earlGrey -g test/test.fasta -s test -o output_dir -t 16`

The archive contains md5 checksums (`checksums.md5`) for all expected output files. Because RepeatModeler uses stochastic sampling, exact output files will vary between runs; users should compare the overall patterns in the high-level count table rather than expecting bitwise-identical results.

Full instructions for running the test and interpreting results were added to the README under a new [Test Your Installation](#test-your-installation) section.

---

### `earlGreyLibConstruct` — quietBar support added

The `-q` flag was fully implemented in `earlGreyLibConstruct` to match the existing behaviour in `earlGrey`:

- **Usage text**: added `-q == Suppress TEstrainer parallel progress bar (yes/no, Default: no, useful for batch/sbatch jobs)`
- **`getopts`**: added `q:` to the option string and `q) quietBar=${OPTARG};;` to the case block
- **`Checks()` validation**: added the same yes/no guard as `earlGrey`, defaulting to `no` with a warning on invalid input
- **`strainer()` and `strainerResume()`**: both TEstrainer call sites now append `$([ "$quietBar" == "yes" ] && echo " -q")`, identical to the `earlGrey` implementation

## Static verification

```bash
cd /data/toby/EarlGrey

# 1. Confirm earlGrey and earlGreyLibConstruct pass ProcNum directly to TEstrainer (no pre-division)
grep -n 'TEstrainer_for_earlGrey.sh.*-t' earlGrey earlGreyLibConstruct

# 2. Confirm strainthreads is no longer present in either script
grep -n 'strainthreads' earlGrey earlGreyLibConstruct  # should return nothing

# 3. Confirm quietBar is fully implemented in earlGreyLibConstruct
grep -n 'quietBar' earlGreyLibConstruct
```

**Expected for check 1:** All TEstrainer invocations in both scripts show `-t ${ProcNum}`.

**Expected for check 2:** No output.

**Expected for check 3:** Six matches — usage text, getopts case, `Checks()` validation (three lines), and both TEstrainer call sites.

---

# v7.2.6 — Fix Dfam 3.9 installation download URL

## Background

The `dfamCheck()` / `Checks()` function in `earlGrey`, `earlGreyLibConstruct`, and `earlGreyAnnotationOnly` generates a helper script (`configure_dfam39.sh`) when it detects that RepeatMasker has not yet been configured with Dfam 3.9 partitions. The generated script contained a `curl` command and a warning message pointing to `https://dfam.org/releases/current/families/FamDB/`. The `current` symlink on the Dfam release server now resolves to a release newer than Dfam 3.9, so users following the generated script would download the wrong database version.

## Changes made

In all three scripts (`earlGrey`, `earlGreyLibConstruct`, `earlGreyAnnotationOnly`), two strings were updated in the `dfamCheck()` / `Checks()` function:

1. The warning message URL changed from:
   ```
   https://dfam.org/releases/current/families/FamDB/
   ```
   to:
   ```
   https://dfam.org/releases/Dfam_3.9/families/FamDB/
   ```

2. The `curl` command written into `configure_dfam39.sh` changed from:
   ```bash
   curl -o 'dfam39_full.#1.h5.gz' 'https://dfam.org/releases/current/families/FamDB/dfam39_full.[0-16].h5.gz'
   ```
   to:
   ```bash
   curl -o 'dfam39_full.#1.h5.gz' 'https://dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.[0-16].h5.gz'
   ```

These are the only changes in this patch. No functional pipeline logic was altered.

## Verification

```bash
# Confirm the old URL is no longer present in any of the three scripts
grep -n 'releases/current' earlGrey earlGreyLibConstruct earlGreyAnnotationOnly  # should return nothing

# Confirm the versioned URL is present in all three scripts
grep -n 'releases/Dfam_3.9' earlGrey earlGreyLibConstruct earlGreyAnnotationOnly
```

**Expected for check 1:** No output.

**Expected for check 2:** Six matches — the warning-message line and the curl-command line in each of the three scripts.

# v7.3.0 update to enable compatibility with Dfam 4.0 and RepeatMasker 4.2.4

First, I need to understand how the new RepeatMasker works with the separate famdb. I also need to understand how this works with RepeatModeler. Then, I need to update the conda packages for RepeatMasker and RepeatModeler, then test a local earl grey version, then update earl grey installation scripts to work with this new format. I also probably need to make a conda recipe for the FamDB tool so that this can be sourced as a dependency in RepeatMasker, RepeatModeler, and Earl Grey. I also need to update the documentation to reflect these changes and ensure that users are aware of the new dependencies when installing Earl Grey.

## Testing FamDB 3.0

Guidance is from: https://github.com/Dfam-consortium/FamDB/blob/master/QUICKSTART.md

```
cd /data/toby/testDIR/
mkdir -p dfam40
cd dfam40

# make a dev environment
mamba create -n earlgrey_730_dev -c conda-forge -c bioconda pip h5py trf rmblast hmmer python=3.10
mamba activate earlgrey_730_dev

# get the latest FamDB release
wget https://github.com/Dfam-consortium/FamDB/archive/refs/tags/3.0.0.tar.gz
tar -zxvf 3.0.0.tar.gz

# this creates FamDB-3.0.0/
# inside here, there is a Libraries directory that I think will be the new place storing the stuff. This could be complex in a conda env as it is inside the env, but maybe we can get it working this way for user friendliness.

# try and download the Dfam components
python /data/toby/testDIR/dfam40/FamDB-3.0.0/utils/download_dfam.py
# this starts an interactive process, which is not great for automation. Let's see if there are options to make this non-interactive. If not, I may need to modify the script or write a wrapper to automate it.
# there doesn't seem to be a non-interactive option in the help section. I will try and see if just specifying the partitions works
python /data/toby/testDIR/dfam40/FamDB-3.0.0/utils/download_dfam.py all
# this doesnt work, let's just try the interactive thing and see what happens
python /data/toby/testDIR/dfam40/FamDB-3.0.0/utils/download_dfam.py

# Okay, so this creates `/data/toby/testDIR/dfam40/FamDB-3.0.0/Libraries/famdb`
# then downloads the gz partitions
# the script is pretty slow, it is downloading the partitions sequentially. I will need to modify this to download in parallel for it to be usable in an automated installation process. For now, I will let it run and see if it completes successfully.
# this downloads the gz files then uncompresses them in `/data/toby/testDIR/dfam40/FamDB-3.0.0/Libraries/famdb/`
# this is at least consistent with 3.9
```

Now, I will see how the new RepeatMasker works with this external format.

```
cd /data/toby/testDIR/dfam40/
wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.2.4.tar.gz
tar -zxvf RepeatMasker-4.2.4.tar.gz
cd RepeatMasker

perl ./configure \
    -trf_prgm /data/toby/miniforge3/envs/earlgrey_730_dev/bin/trf \
    -rmblast_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -hmmer_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -abblast_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -crossmatch_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -default_search_engine rmblast \
    -famdb_dir /data/toby/testDIR/dfam40/FamDB-3.0.0/
```

Now, I will see how the new RepeatModeler works with this external format.
Version 2.0.8 is already on conda, so I will source it and try the config.

```
mamba install repeatmodeler=2.0.8
cd /data/toby/miniforge3/envs/earlgrey_730_dev/share/RepeatModeler/

if [[ "$(uname -s)" == "Linux" ]]; then
    LTR_STRUCTURAL_SEARCH="y"
    CONFIG_OPTIONS+=" \
    -ninja_dir ${PREFIX}/bin"
else
    LTR_STRUCTURAL_SEARCH="n"
    # ninja_dir option not set for osx because package not available in bioconda
fi

# prompt 1: <PRESS ENTER TO CONTINUE>
# prompt 2: confirm path to running perl interpreter
# prompt 3: Configure for LTR structural search [y] or n?
# Answering y for linux; n for osx because NINJA is not available for osx in bioconda
printf "\n\n${LTR_STRUCTURAL_SEARCH}\n" | \
    perl ./configure \
    -mafft_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -ucsctools_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -recon_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -rscout_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -ninja_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -repeatmasker_dir /data/toby/testDIR/dfam40/RepeatMasker/ \
    -trf_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -genometools_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -ltr_retriever_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -rmblast_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -repeatafterme_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/ \
    -cdhit_dir /data/toby/miniforge3/envs/earlgrey_730_dev/bin/
```
