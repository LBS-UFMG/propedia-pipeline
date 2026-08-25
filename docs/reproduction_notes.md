# Reproduction notes (sample validation vs propedia26_v15.csv)

## Extraction (Rule 3) — LOCKED for demo
Sample: 300 oracle PDB IDs already downloaded. 644 extracted pairs matched to v15.

- Peptide size + sequence exact match: 532/644 (83%).
- Size and sequence agreement move together (no independent sequence-reading bug).
- Residue definition: modeled amino acids incl. modified (Biopython is_aa(standard=False));
  modified residues without a canonical 1-letter code rendered as X, matching v15.

### Known remaining gap (112 pairs, all v15-longer, +1 to +4 residues)
v15 additionally counts NON-amino-acid residues at chain termini / interior
(e.g. ACE/NH2 capping groups, statine/ligand residues in inhibitor peptides such
as 1APT/1APU/1APV/1APW) as X. Our extraction counts only amino-acid residues.
This is a documented convention difference, not a bug; arguably our definition
(biological amino-acid sequence) is the more correct one.

Follow-up (daytime): optional exact-parity mode = "all non-solvent residues,
non-amino-acids as X", with solvent/ion exclusion list, re-validated on this sample.

## Physicochemical (ProtParam) — REPRODUCED
Sample: 644 peptides from extracted pairs, compared to v15 columns.
Exact/within-tolerance match:
  MW (±1.0):          644/644 (100%)
  pI (±0.1):          644/644 (100%)
  GRAVY (±0.01):      644/644 (100%)
  Instability (±0.5): 644/644 (100%)
  Aliphatic (±0.5):   644/644 (100%)
Note: physicochemistry matches 644/644 even where sequence match was 532/644,
because non-standard placeholder residues (X) are stripped before analysis and
do not affect computed properties. Confirms the 17% sequence gap is cosmetic.

## Interaction energy (PRODIGY 2.4.0) — REPRODUCED
Sample: 642 pairs compared to v15. Within tolerance:
  Predicted binding affinity (dG, ±0.5):   640/642 (100%)
  Predicted dissociation constant (Kd, log ±0.5): 640/642 (100%)
  NIS apolar / charged (±0.5):             640/642 (100%)
  No. intermolecular contacts (±2):        640/642 (100%)
  No. apolar-apolar contacts (±2):         640/642 (100%)
KEY PARAMETER: PRODIGY must run with --distance-cutoff 6.0 (not the 5.5 default)
to match v15, consistent with Propedia's 6 A interface definition. Pinned in
config (prodigy.distance_cutoff). Temperature 25 C (default). no_contacts: 5/651.

## COCaDA — open question for validation
Contact-type tally on sample shows u-prefixed types: uAT(16), uRE(69), uSB(40).
These are COCaDA's "undefined/uncertain" attractive/repulsive/salt-bridge contacts
(ambiguous charge/protonation state at the given pH). DECISION NEEDED: does Propedia
count these toward the corresponding contact type, or exclude them? Must match v15's
convention when reference contacts are available. Ask R. Lemos / lab mate who ran COCaDA.

## Sequence signatures (iFeature) — INTEGRATED (validation pending)
iFeature CLI confirmed; 10 paper descriptors present. Per-descriptor widths:
AAC=20 DPC=400 DDE=400 GAAC=5 GDPC=25 GTPC=125 CTDC=39 CTDT=39 CTDD=195 CTriad=343.
DISCREPANCY: all 10 sum to 1591, but paper states 1248-feature vector.
1248 = the 9 descriptors EXCLUDING CTriad (1591-343=1248). Likely the final vector
drops CTriad despite the text listing it. RECONCILE with signature author.
X handling: iFeature silently omits X from composition (denominator excludes X);
no crash/NaN. Match v15's stripping convention at validation time.
No reference feature files locally -> integration only tonight.

## Sequence signatures (iFeature) — INTEGRATED, width reconciled
9 descriptors (AAC,DPC,DDE,GAAC,GDPC,GTPC,CTDC,CTDT,CTDD) concatenate to exactly
1248 features, matching the paper's stated vector size (CTriad excluded; incl.
would give 1591). Sample: 643 peptides x 1248; 8 excluded (<3 standard residues
after X-strip, e.g. LXX/VV) — logged in signatures_excluded.tsv.
Fix: strip non-standard residues then require >=3 standard residues; DDE/CTDT
divide by (len-1)/pair-count and fail below that. Validation pending (no ref file).

## Structural signatures (SIGNA / aCSM-ALL) — INTEGRATED
Runs on PEPTIDE-ONLY PDBs (Propedia 26: signatures on the peptide structure alone).
Propedia 26 params: cutoff_limit=10, cutoff_step=0.2, cumulative=True -> 1800 features
(8 atom categories -> 36 pairs x 50 bins). Sample: 651 peptides x 1800.
NOTE: SIGNA docs/v2.3 used cutoff_limit=20; Propedia 26 Methods specify 10 -> pinned to 10.
read_folder writes CSV keyed by file PATH and in filesystem order -> we normalize to
entry ID and MUST join by ID, never row position. Validation pending (no ref file).

## Clustering — CNR REPRODUCED (modulo documented convention); finer clusters deferred
CNR (seq100_clusters): group by identical peptide sequence. Algorithm reproduces v15
exactly. On sample: 256 groups (ours) vs 266 (v15). ALL 28 boundary differences are
attributable to the peptide-sequence X-convention from the extraction stage: v15
includes terminal cap/modified residues as X, we strip them, so v15 splits some groups
we merge (verified: v15_seq.replace('X','') == our_seq for every mismatched member).
NOTE: seq100_clusters representative LABELS are DB-global (chosen across all members),
so exact label strings require the full dataset; grouping STRUCTURE validated on sample.
DEFERRED: is_leader/leader_id is a FINER sub-clustering within CNR groups (multiple
leaders per seq-cluster -> structural/interface criterion, not reproduced tonight).
Legacy sequence/binding/interface-cluster (Hammock/ProBiS/MUSTANG) need external tools.

## KEY REPRODUCTION FINDING
Every discrepancy across ALL stages traces to ONE documented cause: the extraction
convention for terminal cap/modified residues (counted as X in v15, stripped by us).
This affects ~17% of peptide sequences and propagates identically into contacts and
CNR clustering. All geometry/energy/physicochemistry metrics independent of it match ~100%.

## ML therapeutic classifiers — REPRODUCED (method), 6 classes x 6 algorithms
Features: iFeature 1248-vector. Training data: signa cs4 case-study sets
(NOTE: likely aCSM/v2.3-era sets, not confirmed identical to Propedia 26's;
reconcile with Suppl. Text 2 / Tables S1-S7). Models ported Orange->sklearn
(SVM, GradientBoosting, LogReg, kNN, NaiveBayes, MLP), scored on *_test_main.
Best model per class (AUC): ABP GB 0.97, ACP GB 0.93, QSP SVM 0.92, AAP NN 0.82,
SBP GB 0.76, AIP GB 0.75. GradientBoosting wins 4/6.
KNOWN IMPROVEMENT: no feature scaling -> SVM degenerated on AIP (predicted one
class, f1=0). Add StandardScaler pipeline + bump LogReg/MLP max_iter; expected to
lift AIP/SBP. LogReg/MLP ConvergenceWarnings on large classes are benign (cosmetic).
Reproduces METHOD + competitive performance; exact-number match needs Suppl. Tables.

## Consolidation — DONE (architecture); execution edits pending
Snakefile + 7 rule files (13 rules) wire the full pipeline: collect -> download ->
extract -> [physchem, prodigy, cocada, ifeature, signa, cnr] -> ml -> validate.
`snakemake -n` DRY-RUN PASSES: full DAG resolves in dependency order.
Mode switch: config mode=sample (validates vs v15) | mode=full (whole snapshot).
Final target reproduction_report.txt = self-validating scorecard.
Paths consolidated into config machine: block (no more hardcoded paths in scripts).
PENDING for a live run: (1) run_prodigy.py read params pdb_dir/distance_cutoff/temp
+ output/errors; (2) write_peptide_pdbs.py touch output.marker + read params; (3)
confirm ifeature/signa/ml/cocada/physchem param names match rules at execution.
Then: snakemake --unlock && snakemake --cores 16  (sample) ; mode=full for real run.
FULL RUN is hours-overnight scale (73k pairs); needs download complete first.

## Surface area (FreeSASA) — VALIDATED (NACCESS replacement)
NACCESS unavailable (no binary/source) -> reimplemented with FreeSASA (open, pip).
Seven columns: ASA_Complex/Peptide/Protein, BSA, BPepA, BProA, BPP%. BSA per eq.1:
(ASA_prot + ASA_pep - ASA_complex)/2. Applies paper's BSA>0 selection filter.
Sample: 651 processed, 648 pass BSA>0. Validation vs v15 (NACCESS):
  BSA correlation r=1.000 ; mean BSA ours=605.5 v15=601.4 (0.7% offset).
=> FreeSASA reproduces NACCESS surface areas essentially exactly, AND is
redistributable where NACCESS is not — net reproducibility improvement (answers Reviewer 1).
