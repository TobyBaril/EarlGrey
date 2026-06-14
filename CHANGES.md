# Changes

Summary of additions and changes on this branch (`hyphaltip_explore_claudeImprove`)
compared to `main`.

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
