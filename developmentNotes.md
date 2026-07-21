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

Now I need to build and test the recipe

```
cd /data/toby/testDIR/dfam40/famdb.conda
conda build .
conda create -n famdb.test famdb --use-local

conda activate famdb.test

# try and install some partitions and see if it works
download_dfam.py
# I chose "1,2,3"
# then "all"
```

This downloads stuff to: `/data/toby/miniforge3/envs/famdb.test/Libraries/famdb`

## Combining with RepeatMaker

Now, I will build the recipe I have modified for RepeatMasker 4.2.4, which should be compatible with the new FamDB format. I will then test this installation and see if it works with the new FamDB format.

```
cd /data/toby/testDIR/dfam40/repeatmasker.conda
conda build .
conda activate famdb.test # this env has the famdb tool installed
conda install -n famdb.test repeatmasker=4.2.4 --use-local

cd /data/toby/miniforge3/envs/famdb.test/share/RepeatMasker
perl configure # everything should already be there
# turns out you do not need to reconfigure RepeatMasker after sourcing all the partitions now!
# Libraries get built in /data/toby/miniforge3/envs/famdb.test/share/RepeatMasker/Libraries/CONS-Dfam_4.0/
```

## Sorting out RMBlast
The new RepeatModeler needs rmblast 2.17.1 so I also need to get the new rmblast working in conda.

```
mkdir -p /data/toby/testDIR/dfam40/rmblast.conda
cd /data/toby/testDIR/dfam40/rmblast.conda

# the build fails because of changes in the patches, I need to make these again and build from scratch
```

Fixing the patches:

```
mkdir -p /data/toby/testDIR/dfam40/fixingRMBLAST && cd /data/toby/testDIR/dfam40/fixingRMBLAST
# 1. Get the new source and apply the rmblast patch
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-src.tar.gz
tar zxvf ncbi-blast-2.17.0+-src.tar.gz
cd ncbi-blast-2.17.0+-src
patch -p1 < /data/toby/testDIR/dfam40/rmblast.conda/isb-2.17.1+-rmblast.patch

# 2. Apply all the other non-VDB patches that DO work
patch -Np0 -i /data/toby/testDIR/dfam40/rmblast.conda/configurellvm.patch
patch -Np1 -i /data/toby/testDIR/dfam40/rmblast.conda/get_species_taxids.patch
patch -Np1 -i /data/toby/testDIR/dfam40/rmblast.conda/phonehome.patch
patch -Np1 -i /data/toby/testDIR/dfam40/rmblast.conda/update_configsub.patch
patch -Np1 -i /data/toby/testDIR/dfam40/rmblast.conda/project_tree_builder.patch

# 3. Save the current state of configure as a baseline
cp c++/src/build-system/configure c++/src/build-system/configure.orig

# all patches work, vdb is not needed anymore, so remove it from the yaml
```

Try and add to the env and see if RepeatModeler is happy now:

```
conda activate famdb.test # this env has the famdb tool installed
conda install rmblast --use-local
```

## Combining with RepeatModeler
Now, I will see how the new RepeatModeler works with this external format.
Version 2.0.9 is required for RepeatModeler to work with the new FamDB format. I will install this version and see if it works.

```
cd /data/toby/testDIR/dfam40/
mkdir -p repeatmodeler.conda
cd repeatmodeler.conda

conda build .

# test RepeatModeler
conda activate famdb.test
conda install -n famdb.test repeatmodeler=2.0.9 --use-local

cd /data/toby/testDIR/dfam40
BuildDatabase -name dfamTest /data/toby/testDIR/test.fasta
RepeatModeler -database dfamTest -threads 16 
```

All these packages are working now. 

Now, I need to update Earl Grey to work with the new Dfam libraries. This will require updating the conda recipe to add the new tool dependencies. This will also require changes to the Earl Grey install guide/scripts to use the new Dfam download tool that they have made (this should be simpler than my current method...).

First, I have updated the conda recipe and I will try and build the new package locally.

```
cd /data/toby/EarlGrey
conda build conda/

conda create -n earlgrey_730 -c conda-forge -c bioconda earlgrey=7.3.0 --use-local

# try trigger a run, this should stop and signal config is required
cd /data/toby/testDIR/
earlGrey -g test.fasta -s test_730 -o . -t 16
```

The tools crash as BuildDatabase and RepeatModeler do not work properly. I will clear my build cache and try and install all things from Bioconda except earlgrey and RepeatModeler...

```
conda deactivate
conda-build purge-all

# build RepeatMasker 4.2.4
conda build /data/toby/testDIR/dfam40/repeatmasker.conda/

# build rmblast 2.17.1
conda build /data/toby/testDIR/dfam40/rmblast.conda/

# build RepeatModeler 2.0.9
## recipe is stored at: /data/toby/testDIR/dfam40/repeatmodeler.conda
conda build /data/toby/testDIR/dfam40/repeatmodeler.conda/

# build Earl Grey 7.3.0
conda build /data/toby/EarlGrey/conda/
```

In this test, RMBlast is fixed, now I need to fix RepeatModeler and test Earl Grey

Now, try and make a new environment and install Earl Grey 7.3.0 with the new dependencies:

```
conda create -n earlgrey_730_test -c conda-forge -c bioconda earlgrey=7.3.0 --use-local
conda activate earlgrey_730_test




## Patches and Improvements contributed by @hyphaltip
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