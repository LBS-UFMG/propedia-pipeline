# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A **Snakemake** pipeline that rebuilds the **Propedia 26** protein–peptide interaction database
from a PDB snapshot. Two datasets are produced: **pep-pro** (`propedia.csv`, one peptide vs one
protein chain) and **multipro** (`multipro_final.csv`, one peptide vs ≥2 protein chains). The
pipeline replaces an older ad-hoc script-based build; it is the canonical way to (re)build Propedia
and is the subject of a NAR Database revision. `README.md` is the user-facing guide and
`docs/reproduction_notes.md` records every design decision and validation result — read those for
depth; this file is the orientation + the non-obvious gotchas.

Linux only (case-sensitive chain IDs collide on NTFS). Run everything from the repo root.

## Commands

```bash
snakemake -n --cores 1                 # DRY RUN — always do this after editing rules
snakemake --cores <N>                  # sample mode (default) — builds + self-validates
snakemake --cores <N> --config mode=full --rerun-incomplete --keep-going   # production
snakemake --cores <N> package          # build the website file tree (LOCAL only, no deploy)
snakemake --cores <N> --config mode=full site    # ONE-COMMAND production build: dataset +
                                       #   packaged web tree (data/db/…) + bulk zips, ready to rsync
scripts/deploy_site.sh user@host:/…/public/data  # PURE TRANSFER of the `site` tree (dest never
                                       #   hardcoded; also honors $PROPEDIA_DEPLOY_DEST; -n = dry run)
snakemake --cores 8 --config mode=full release   # dated immutable release snapshot
snakemake --unlock                     # after any hard-killed run, before resuming
python -m py_compile workflow/scripts/<x>.py      # there is no unit-test suite; this + the
python tests/smoke_test_pisa.py --ccp4-dir <ccp4> # PISA smoke test are the checks
```

- `--cores <N>` is the whole parallelism knob together with `compute.threads` in the config;
  set `--cores` to your allowed core budget. The heavy per-entry stages honor `compute.threads`.
- Shrinking the sample: `--config sample='{size: 30, seed: 42}'`. This only re-runs `build_sample`;
  if `pairs.tsv` is already cached from a bigger sample, Snakemake reuses it and the override does
  **not** shrink the data. For a genuinely small run: `rm -rf results/sample state/sample` first.
- The final target is `results/<mode>/reproduction_report.txt` — it diffs each stage against the
  reference (v15) oracle. `package`, `site`, and `release` are separate explicit targets (`site`
  = report + `package` + `zips`; kept out of `rule all` so sample runs don't build the multi-GB zips).
- **Deploy = pure transfer.** `package.web_mode_name: "db"` writes the tree straight into `data/db/`
  (the site's `$modo`), so `snakemake site` then `scripts/deploy_site.sh <dest>` needs no rename step.
  The packager also emits `data/pdb/total_contacts2.txt` (home-page counters, incl. the seq-dedup
  unique count) so the site's stats stay in sync every build.

## Architecture

**DAG:** `workflow/Snakefile` includes `workflow/rules/01…11_*.smk`; each rule runs a script in
`workflow/scripts/`. Flow: collect PDB IDs → download CIFs → `extract_pairs` (per-pair mmCIF via
Biopython `MMCIFIO`) → feature fan-out (physchem, freesasa, prodigy, cocada, ifeature, signa,
interface, metadata, pisa) → clusters + ML/therapeutic → `assemble` → multipro (Phase 1 grouping →
Phase 2 recompute → `multipro_assemble`) → `validate` → `package`/`release`.

**Three output trees, and the distinction is load-bearing:**
- `state/<mode>/` — **DURABLE** database memory: per-entry checkpoints + extracted structures
  (`cif/`, `pep_pdb/`, `multipro_cif/`, `cocada/`). Preserve & back up; this is what makes updates
  incremental. Git-ignored.
- `results/<mode>/` — **DERIVED** TSVs/CSVs + report. Disposable, rebuilt from `state/`.
- `releases/propedia-<date>/` — immutable dated snapshots + manifest.

**Per-entry checkpointing (`workflow/scripts/checkpoint.py`)** — the heavy stages (extract,
prodigy, cocada, freesasa, interface, metadata, signa, pisa, multipro_*) define a picklable
`worker(item)` and run a bounded process pool. One JSON per entry is persisted to a work dir
Snakemake does NOT track, so a killed run resumes per-entry. Key rules:
- The work dir is namespaced by a hash of the stage's `VERSION` + params. **Bump a script's
  `VERSION` to force a clean recompute of that stage** (or change its params); old checkpoints are
  then ignored.
- A status starting `retry:` is transient — NOT persisted, retried next run. Any other status
  (`ok`, `no_output`, `timeout`, …) IS persisted (never recomputed). Use `retry:` only for truly
  transient failures; use a persisted status for deterministic per-entry outcomes.

**Config:** `config/config.yaml`. Only the top `machine:` block is machine-specific (external tool
paths, CIF dir, oracle CSVs). Everything else is pipeline science — don't change without intent.
`mode: sample|full`.

## Non-obvious conventions & gotchas (learned the hard way)

- **mmCIF intermediate.** Extraction writes ATOM-only mmCIF via Biopython `MMCIFIO` (lossless for
  >99,999 atoms, multi-char chain IDs, wide numbering). Downstream tools key on `auth_asym_id`.
- **COCaDA cannot parse `MMCIFIO`-generated CIFs** (`'NoneType' object is not iterable`). Both
  `run_cocada.py` (pep-pro) and `multipro_cocada.py` therefore run COCaDA on the **original
  downloaded RCSB CIF**, not the extracted one. COCaDA selects on **`label_asym_id`**, so chain IDs
  are translated `auth_asym_id → label_asym_id` (see `auth_to_polymer_label`) before `-c`.
- **SIGNA calls bare `exit()`** on a degenerate peptide ("No match found.") — which would kill the
  host process mid-batch. `run_signa.py` calls SIGNA's `read()` **per file** and catches
  `SystemExit`, skipping only the bad peptide. Do not switch back to `read_folder`.
- **PISA is optional.** `pisa` stage needs CCP4 (`machine.ccp4_dir`). If absent it honors
  `pisa.on_missing`: `skip` (default; blank cols, `pisa_status=pisa_unavailable`) or `error` (use
  for production). X-ray only (cryo-EM/NMR → `not_applicable`). The derived `pisa_interface_class`
  (biological / crystal-packing) comes from PISA's CSS vs `pisa.css_biological_threshold` (0.5),
  applied at MERGE time so re-labeling never re-runs PISA. This annotation answers the reviewer
  concern about crystal-packing artifacts.
- **X-residue convention (v15-style, restored).** Sequences include every **polymer** residue —
  terminal caps and modified residues become `X`, matching v15. `extract_pairs.polymer_set()`
  reads polymer membership from `_atom_site.label_seq_id` (numeric = polymer; `.`/`?` =
  water/ligand/ion, excluded); `modeled_aa` keeps those (falling back to the amino-acid test only
  when a cif has no `label_seq_id`). Earlier versions dropped caps (`is_aa`), which caused a
  documented ~17% peptide-sequence gap and ~69% report match; that is now ~100%. Note: only the
  SEQUENCE path carries the caps — the extracted per-pair cif still writes `is_aa` residues
  (`PairSelect`), so `PEPTIDE_SIZE`/seq (padded) do not equal the cif's residue count (harmless;
  nothing indexes seq against cif). The size filter (`pep 2–50`) applies to the padded length.
- **`--config` values arrive as strings** when overriding nested keys; cast (e.g. `int(...)`) —
  see `build_sample.py`.
- **Structure size gate (`cutoffs.max_atoms_per_structure`)** is a MEMBERSHIP filter, not a
  compute/checkpoint param. `extract_pairs` cheaply counts `_atom_site` ATOM/HETATM rows
  (no Biopython parse), caches them in `state/<mode>/atom_counts.tsv`, and drops over-limit
  PDBs from `pairs.tsv` (so they never parse or cascade). It is deliberately NOT in the
  extract checkpoint signature: already-computed entries survive a limit change, and RAISING
  the limit re-admits previously-skipped structures as new work without recompute. Do NOT
  add it to the `checkpoint.namespace` params, and do NOT bump extract `VERSION` for it.
- **Legacy clusters are inherited, not recomputed** (Hammock/MUSTANG/ProBiS). New entries get blank
  `sequence-/interface-/binding-cluster` + `is_leader`/`leader_id` until those tools are re-run.
  The `seq100`/CNR cluster IS recomputed.
- **Column order for the website** is exact and position-based. Pep-pro per-entry files use the 93-col
  v17 order in `config/columns_peppro_v17.txt`; multipro is already emitted in the correct v4 order.
  `make_package.py` matches target→source columns case-insensitively (site `PISA_*` ↔ our `pisa_*`).
  Packaging writes a LOCAL `web_dir` tree only — it never deploys.

## External tools (paths in `config/config.yaml` `machine:`)

`prodigy-prot` (via PATH), COCaDA, SIGNA, iFeature (cloned repos), CCP4≥9 for PISA (optional).
Activate the Python env before running so subprocess tools resolve on `PATH`.
