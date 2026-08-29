# Propedia Rebuild Pipeline

A reproducible [Snakemake](https://snakemake.github.io/) pipeline that rebuilds the
Propedia protein–peptide interaction database from a Protein Data Bank (PDB) snapshot.
It collects structures, extracts protein–peptide pairs, computes every feature family
(physicochemistry, surface area, interaction energy, contacts, sequence & structural
signatures), clusters, trains therapeutic-peptide classifiers, and validates the result
against a reference release.

> **Status (read before running):** the pipeline architecture is complete and its DAG
> is verified (`snakemake -n` passes). A few per-script parameter wirings still need to
> be finished before a clean one-command run — see **[Current state](#current-state)**.
> This is a working research pipeline, not yet a turnkey release.

---

## 1. Requirements

- **OS:** Linux (developed on Ubuntu/WSL2). Use a **native Linux filesystem**, not a
  Windows-mounted drive — PDB/chain IDs are case-sensitive and collide on NTFS.
- **Python** ≥ 3.10, with:
  ```
  pip install snakemake biopython freesasa prodigy-prot scikit-learn scipy numpy
  ```
- **External tools** (clone these; paths go in the config):
  - COCaDA — `git clone https://github.com/LBS-UFMG/COCaDA` (checkout the tag you want; ≥ v1.5 for the `-c/-inter/-ph/-s` flags)
  - SIGNA — `git clone https://github.com/LBS-UFMG/signa`
  - iFeature — `git clone https://github.com/Superzchen/iFeature`
  - **CCP4 ≥ 9** (for PISA) — provides `bin/pisa` + `share/pisa/pisa.cfg`; set `machine.ccp4_dir`.
    PISA powers the **biological-vs-crystal-packing** interface annotation
    (`pisa_interface_class`). It is **optional**: if CCP4 is absent, the `pisa` stage
    honours `pisa.on_missing` — `"skip"` (default) builds without PISA (blank PISA
    columns, `pisa_status='pisa_unavailable'`); set it to `"error"` for a production
    release to require PISA. PISA runs on **X-ray** entries only (no crystal lattice for
    cryo-EM/NMR → `pisa_interface_class='not_applicable'`).
- **Disk:** ~30–60 GB for the CIF download (full mode).
- **Reference DB** (for validation): a Propedia release CSV (e.g. `propedia26_v15.csv`),
  semicolon-delimited.

---

## 2. Directory layout

```
.
├── config/
│   └── config.yaml            # ALL settings; edit the `machine:` block for your paths
├── workflow/
│   ├── Snakefile              # pipeline entry point (includes the rule files)
│   ├── rules/                 # 01_collect … 07_validate  (the DAG)
│   └── scripts/               # the Python each rule runs
├── docs/
│   └── reproduction_notes.md  # detailed technical log of every stage + known caveats
├── state/                     # DURABLE db memory: checkpoints + extracted structures
│   └── <mode>/                #   (git-ignored; PRESERVE & back up — enables incremental updates)
├── results/                   # DERIVED outputs (TSVs/CSVs/report); disposable, rebuilt from state
│   └── <mode>/
└── releases/                  # dated release snapshots + provenance manifests (git-ignored)
    └── propedia-<date>/
```

`config/` and `workflow/` are the code (version them). `results/` is disposable (rebuilt
from `state/`). `state/` is the durable database memory — preserve it between releases so
an update only fills in new PDB entries. A release is fully described by its manifest
(config + git commit + snapshot date).

---

## 3. Configure (the only machine-specific step)

Open `config/config.yaml` and edit the top **`machine:`** block to your paths:

```yaml
machine:
  cif_dir:      "/path/with/space/for/cifs"        # CIF download target (~50 GB)
  oracle_csv:   "/path/to/propedia26_v15.csv"      # reference DB for validation
  cocada_dir:   "~/COCaDA"                          # cloned COCaDA repo
  signa_dir:    "~/signa/Signa"                     # dir containing signa.py
  ifeature_dir: "~/iFeature"                        # cloned iFeature repo
  ml_train_dir: "~/signa/docs/case_studies/cs4/csv" # therapeutic-class train/test sets
```

Everything below `machine:` is pipeline logic (cutoffs, tool parameters) — leave it
unless you intend to change the science. Key pinned parameters, and why:

| Setting | Value | Why |
|---|---|---|
| `prodigy.distance_cutoff` | 6.0 | matches Propedia's 6 Å interface (not PRODIGY's 5.5 default) |
| `signatures.acsm_cutoff_limit` | 10 | Propedia 26 aCSM setting → 1800-feature vector |
| iFeature descriptors | 9 (no CTriad) | → the paper's 1248-feature vector |
| `cutoffs.*` | — | distance/BSA/length selection criteria |

---

## 4. Run

Pick a mode with `mode:` in the config (or override on the command line).

**Sample mode** (fast; rebuilds a sample and validates it against the reference DB):
```bash
snakemake --cores 16                      # mode defaults to "sample"
```

**Full mode** (rebuilds the whole snapshot — this is the production run):
```bash
snakemake --cores 16 --config mode=full
```

Always dry-run first to see the plan without executing anything:
```bash
snakemake -n --cores 1
```

The final target is a self-validating report:
```
results/<mode>/reproduction_report.txt
```
which diffs each stage against the reference DB and prints the reproduction scorecard.

### Long runs: parallelism & resuming (read before a full run)

The heavy per-entry stages (extract, prodigy, cocada, freesasa, interface, metadata,
signa, ifeature, multipro) **checkpoint every entry** and run a bounded process pool.
Two consequences:

- **Parallelism** is one knob: `compute.threads` in the config (default 8). Set it to
  your **physical** core budget and run `snakemake --cores <logical threads>`; snakemake
  governs total concurrency so co-scheduled stages don't oversubscribe. Setting threads
  to the physical core count (not the SMT thread count) keeps each CPU-bound stage at its
  sweet spot, while `--cores` at the thread count lets snakemake co-schedule two stages
  (e.g. CPU-bound freesasa ‖ wait-heavy cocada) to fill the SMT siblings. Example, on a
  box where you may use **32 cores / 64 threads**:
  ```bash
  # config: compute.threads: 32
  tmux new -s propedia            # so an SSH drop doesn't kill the run
  snakemake --cores 64 --config mode=full --rerun-incomplete --keep-going
  ```
  (Budget RAM for the pool — see below. Here ~64 workers × up to ~2 GB stays well under
  a 384 GB allowance.)
- **Resumable at per-entry granularity.** If the run is interrupted (kill, crash, SSH
  drop, power), just re-run the **same command**. Finished stages are skipped; an
  interrupted stage resumes from its checkpoints and only computes the entries it hadn't
  finished — it does **not** recompute from zero. This state lives in `state/<mode>/`
  (see **Durable state** below).
  - After a hard kill, run `snakemake --unlock` once before resuming.
  - To force a clean recompute of a stage, delete `state/<mode>/checkpoint/<stage>/`
    (or bump the stage script's `VERSION`, which auto-invalidates its checkpoints).
    Changing a stage's params also auto-invalidates.

### Durable state, updates, and releases

Two output trees, and the distinction matters:

- **`state/<mode>/`** — the **durable database memory**: per-entry checkpoints plus the
  extracted per-pair structures (`cif/`, `pep_pdb/`, `multipro_cif/`, `cocada/`).
  **Preserve and back this up between releases.** It is what lets an update process only
  the *new* PDB entries. (git-ignored; not disposable.)
- **`results/<mode>/`** — **derived** TSVs/CSVs and the report. Cheap to rebuild from
  `state/`; safe to wipe.

**Updating the database** (the production workflow): bump `pdb_snapshot_date` in the
config and re-run full mode with the *same code*. The download fetches only newly
released CIFs; every stage finds existing entries in `state/` and **skips them**,
computing only the new PDB entries; the CSVs are rebuilt with old + new together. Keep
`state/full/` intact and this is incremental — wipe it and it rebuilds from scratch
(still correct, just slow).

> Changing a stage's *logic* (bumping its `VERSION`) intentionally recomputes that stage
> for **all** entries, so a release is internally consistent rather than a mix of old and
> new logic.
>
> **Known gap:** the legacy sequence/interface/binding clusters (`is_leader`, `leader_id`)
> are inherited from frozen tool outputs and are **not** recomputed, so brand-new entries
> get blank legacy-cluster columns until those tools are re-run. (The `seq100`/CNR cluster
> *is* recomputed.)

**Cutting a release:** after a run finishes, snapshot the deliverables into a dated,
immutable directory with a provenance manifest (snapshot date, git commit, tool config,
entry counts):
```bash
snakemake --cores 8 --config mode=full release
```
→ `releases/propedia-<snapshot_date>/{propedia.csv, multipro_final.csv, reproduction_report.txt, manifest.json}`

**Other notes:**
- The **CIF download** is separately incremental (skips validated files); it too resumes.
- **Full mode is a long run** (hours to overnight). Give it enough RAM for the pool: the
  parse-heavy stages hold one structure per worker (up to ~2 GB for the largest modern
  complexes), so budget roughly `threads × 2 GB` headroom.
- Subprocess tools (`prodigy`, COCaDA's `python3`, iFeature) resolve via `PATH`, so
  **activate the environment first** (`source .venv/bin/activate`) or run snakemake with
  the env on `PATH`.

---

## 5. What each stage does

| Rule | Output | Tool |
|---|---|---|
| `fetch_candidate_ids` | candidate PDB IDs | RCSB Search API |
| `download_cifs` | gzipped mmCIF files | RCSB file service |
| `extract_pairs` | `pairs.tsv` (+ per-pair mmCIF) | Biopython |
| `physchem` | physicochemical properties | Biopython ProtParam |
| `freesasa` | ASA/BSA surface areas | FreeSASA |
| `prodigy` | binding ΔG, Kd, contacts | PRODIGY |
| `cocada` | interatomic contacts | COCaDA |
| `ifeature` | sequence signatures (1248) | iFeature |
| `peptide_pdbs` → `signa` | structural signatures (1800) | SIGNA (aCSM) |
| `cnr_cluster` | 100%-identity peptide clusters | — |
| `ml_classifiers` | 6 therapeutic-class models | scikit-learn |
| `validate` | reproduction scorecard | — |

See `docs/reproduction_notes.md` for the validation results and every design decision.

---

## 6. Current state

**Verified:** the full DAG resolves (`snakemake -n` passes). On a validation sample,
extraction, physicochemistry, surface area (FreeSASA vs NACCESS, r≈1.0), and PRODIGY
reproduce the reference release; signatures, contacts, clustering, and ML are integrated
with correct outputs. Details and caveats in `docs/reproduction_notes.md`.

**Unfinished — needs doing before a clean full run:**
1. **Wire two script params to the rules:** `run_prodigy.py` (read `pdb_dir`,
   `distance_cutoff`, `temperature` from `snakemake.params`) and `write_peptide_pdbs.py`
   (touch its `output.marker`). Others were run via a test harness and should be confirmed
   against their rule signatures at first execution.
2. **Not yet built (paper scope not fully covered):**
   - **Multipro dataset** (peptide vs. ≥2 protein chains) — a whole second dataset.
   - **Interface-residue** annotation column.
   - A few **physicochemical columns** (atomic formula, atom count, extinction coefficient).
   - Some **PDB metadata** columns (resolution, method, deposition date, title, organism).
   - **Legacy clusters** (Hammock / MUSTANG / ProBiS) and the **PepBDB/PepBind** comparison.
3. **Cleanup:** pin tool versions, add feature scaling to the ML step, vendor COCaDA at a
   fixed tag instead of a config path.

**Known caveat that explains most discrepancies:** the reference DB counts terminal
cap / modified residues as `X`; this pipeline counts only standard amino acids. This
affects ~17% of peptide sequences and propagates identically into contacts and clustering.
All metrics independent of it match ~100%. This is a deliberate, documented convention.

---

## 7. Reproducing / updating the database

To update Propedia with newer structures: set `pdb_snapshot_date` in the config to the new
date, set `mode: full`, and run. The pipeline pulls everything released up to that date,
processes only what's new (the download and manifest are incremental), and rebuilds.
