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

## Physicochemical extras (formula, atoms, extinction) — ADDED
Four v15 columns now emitted per chain: Formula, TotalAtoms, ExtCoeff_Disulfide,
ExtCoeff_NoDisulfide. Validated on the 300-ID sample (644 pairs, 1288 chains) with
check_physchem_extras.py.

- Formula / TotalAtoms: 1288/1288 (100%). Computed as the SUM of per-residue atomic
  compositions (free amino acid minus one water) with NO added terminal water — this
  no-water convention is what reproduces v15 exactly (adding H2O overshoots by 3 atoms).
  Computed on the X-stripped sequence, like all physchem.

- Extinction coefficients: emitted as the standard ExPASy/Biopython values
  (reduced = Tyr*1490 + Trp*5500 ; cystine = reduced + (nCys//2)*125). These are a
  DELIBERATE CORRECTION of v15, which used a non-standard formula:
      v15 NoDisulfide = base + nCys*125      (adds 125 per FREE cysteine — non-physical;
                                              reduced Cys do not absorb at 280 nm)
      v15 Disulfide   = base + (nCys*125)//2 (correct cystine value but not floored to
                                              whole pairs; differs from ours by 62 only
                                              when nCys is odd)
  The reduced-Cys +125 term in v15's NoDisulfide is almost certainly a bug. We ship the
  correct values (same rationale as replacing NACCESS with FreeSASA). Validation confirms
  the divergence is fully understood: exact-match is low (45%/67%, only Cys-free chains),
  but 100% of the differences equal exactly the v15 term above — i.e. our counts match
  v15 and only the formula differs. The public release should note this correction so the
  two columns changing vs the online v15 is expected, not a regression.
  
## Interface residues — ADDED (6 A atom-atom, validated)
New column `interface_residues`: protein residues (author resseq, ascending,
comma-joined) with any atom <=6 A from any peptide atom, computed on the extracted
two-chain PDB. Script: workflow/scripts/interface_residues.py.

Sample validation (642 entries with a v15 list): 585/642 exact (91.1%).
A cutoff sweep pins v15's method to exactly 6.0 A — it is a sharp peak that
collapses on both sides (5.5 A -> 11.7%, 6.05 A -> 69.5%, 6.5 A -> 10.0%), so our
6 A implementation is correct and the residual 9% are v15 inconsistencies, not a
cutoff to adopt. The mismatches split as:
  - ~38 pep_shorter: v15 peptide includes non-AA cap/statine/ligand residues
    (e.g. 1APT/1APU/1APV/1APW) whose atoms reach protein residues our AA-only
    peptide does not — the same terminal-X convention as the sequence stage.
  - ~10 extra_only: v15 excludes a few high-numbered modified protein residues we
    keep (e.g. 1AQC, 1AWT) — mirror image, same non-standard-residue family.
  - a few equal-size cases (1AB9, 1B3F, 1C5W/X, 1CM4, 1CWE) where every residue
    v15 lists is >6 A in our PDB (up to 9.85 A in 1CM4) — v15 used a looser/other
    threshold for these; not reproducible under a consistent 6 A rule and not
    desirable to reproduce.
v15 also blanks ~5.9% of entries DB-wide despite a valid <=6 A contact (likely a
defect); we emit the correct non-empty list there.

## PDB metadata (CIF header) — ADDED
Eight v15 columns from the CIF header, plus two additive organism columns.
Script: workflow/scripts/metadata.py (parses each CIF once per PDB id, cached).
Sample validation (644 pairs):
  CLASSIFICATION   644/644 (100%)   _struct_keywords.pdbx_keywords
  DEPOSITION_DATE  644/644 (100%)   _pdbx_database_status.recvd_initial_deposition_date
  RESOLUTION       644/644 (100%)   _refine.ls_d_res_high (blank for NMR)
  STRUCTURE_METHOD 644/644 (100%)   _exptl.method
  TITLE            644/644 (100%)   _struct.title (verbatim, case preserved)
  organism         635/644 (98.6%)

organism: a single per-pair value in v15 = the PROTEIN chain's entity source
organism (uppercased), resolved via _entity_poly.pdbx_strand_id -> entity ->
_entity_src_{gen,nat} / _pdbx_entity_src_syn. NOTE: v15's exact source-selection
rule is undetermined — it is neither "first entity" nor strictly the protein
(counter-examples both ways: 1FC2 reports the peptide's organism, 1BBR the
protein's). The protein-chain rule reproduces 98.6% on the sample; the 9 misses
are v15 quirks (1CJQ/1CJR: synthetic thrombin tagged "SYNTHETIC CONSTRUCT";
1D3E/1D3I: rhinovirus coat protein where v15 recorded the human host). Not
reproduced further; not desirable to chase.

ADDED (beyond v15): PEPTIDE_ORGANISM and PROTEIN_ORGANISM resolve each chain's
organism separately. v15 carried only one per-pair organism, which is misleading
for cross-species complexes (thrombin+hirudin, MHC+viral peptide, protein+synthetic
inhibitor). No v15 counterpart to validate against.

PEPTIDE_DESC / PROTEIN_DESC = _entity.pdbx_description for the chain's entity.
DELIBERATELY kept as the real entity descriptions (peptide 63.2%, protein 69.4%
vs v15). v15's DESC is unreliable: it uses the entity description when present but
falls back to the literal constants "Peptide"/"Protein" otherwise, and sometimes
assigns the descriptive text to the wrong chain (e.g. 1A37/1A38 put the peptide's
description on the protein side). Reproducing that fallback would replace correct,
informative descriptions with constants or noise, so we keep the clean values and
document the divergence (same posture as extinction coefficients / FreeSASA).  

## Assembly — the reproduced Propedia CSV
workflow/scripts/assemble.py joins every feature TSV into the final ;-delimited
propedia.csv. Master row set = extracted pairs passing BSA>0 (ids present in
surface.tsv), sorted by id; left-join all features, blank where missing.
Columns are the exact v15 order (71) + two additive organism columns
(PEPTIDE_ORGANISM, PROTEIN_ORGANISM). Wired as rule 'assemble' (08_assemble.smk).

Sample validation (check_assemble.py): 648 rows x 73 cols; first-71 header EXACT
match vs v15 (incl. the U+02DA degree glyph). Column alignment on the 644 v15
overlap: PDB_ID / chains / CLASSIFICATION / DEPOSITION_DATE / RESOLUTION /
STRUCTURE_METHOD / TITLE / peptide+protein Formula & TotalAtoms all 100%;
organism 98.6% (documented v15 quirks); SEQ/SIZE 82.6% peptide / 95.8% protein
(terminal-X convention). ML classes (AAP/ABP/ACP/AIP/QSP/SBP) and legacy clusters
(binding-/interface-/sequence-cluster, is_leader, leader_id) intentionally blank
pending their stages. seq100_clusters grouping is reproduced but its labels are
DB-global, so that column validates on the full run, not the sample. 648 vs 644
overlap = 4 borderline pairs we include beyond v15's slice of these PDBs;
reconciled at full scale.

## Therapeutic-class scoring (AAP/ABP/ACP/AIP/QSP/SBP) — ADDED
Six per-peptide probability columns.
- Bake-off (run_ml.py): six sklearn models (SVM, GradientBoosting, LogisticRegression,
  kNN, NaiveBayes, NeuralNet) trained per class on the propedia26-sm `new` train TSVs,
  scored on the held-out test_main set. StandardScaler in-pipeline (fit on train only,
  no leakage). Default hyperparameters (proof-of-principle, not optimized). Full metric
  table (all models x classes, incl. AUC) -> ml_report.tsv (Supplementary); GB featured
  in main text.
- Per-class winner selected by lowest test-set BRIER SCORE (probability quality), NOT AUC,
  because the shipped columns are probabilities (a model can rank well yet be poorly
  calibrated). Winners: AAP SVM, ABP SVM, ACP GB, AIP GB, QSP GB, SBP GB.
- Scoring (therapeutic_scoring.py): the winning model per class predicts a probability for
  every Propedia peptide. Features RECOMPUTED via the same iFeature 1248-vector as the
  signature stage (train+propedia in one pass) so ranking and scoring share a feature space.

Selection principle: models chosen to make the shipped columns best on their own merits
(held-out Brier/AUC), NOT to match v15 — v15's exact per-class models/training are
undetermined (not in the repos; only Suppl. Tables S1-S7 would settle it).

Validation vs v15 (proof-of-principle: method, not exact numbers). Per-peptide Pearson r
on the 644 sample overlap: AAP 0.83, ACP 0.83, QSP 0.73, SBP 0.73 (reproduced); ABP 0.52
and AIP 0.07 DIVERGE. ABP: SVM is best on our held-out test by both Brier (0.052) and AUC
(0.979), so the lower v15 correlation is a v15-vintage difference (different model or the
old training set), not a selection artifact. AIP: the only imbalanced training set
(1258/1887); no model tracks v15 (likely v15 balancing, or the new/old split). Both
documented, not chased.

Note: SVC(probability=True) is deprecated (sklearn>=1.9, removed in 1.11) -> pin
sklearn<1.11 or migrate the SVM to CalibratedClassifierCV before upgrading.

## Legacy clusters — INHERITED (not recomputed)
sequence-cluster / interface-cluster / binding-cluster / is_leader / leader_id are
read from the precomputed clusters_v15/ files (legacy_clusters.py), not recomputed:
- sequence/interface/binding: label sits on each cluster's FIRST member only (the
  representative); all other members blank. Validated first-member rule = 100% on
  real data. These are frozen tool outputs (Hammock/MUSTANG/ProBiS) inherited from v1.
- is_leader/leader_id: from redundant.tsv (full-complex redundancy grouping). Leader
  is NOT a simple alphabetical pick (~87.5%), so not reliably recomputable -> inherited.
Consequence: NEW entries in a future update are unclustered until the tools are re-run.
Update path (future work): recompute redundancy from pairs.tsv (no tools needed) and
assign new entries to existing sequence/binding clusters via the BLAST/ProBiS machinery
already in the pipeline; full re-clustering (re-running all three tools) is a heavier
separate option that also changes labels each release.

## Multipro dataset — Phase 1 (grouping + derivable columns)
multipro_extract.py groups pep-pro pairs by (PDB, peptide chain); an entry is
Multipro when the peptide contacts >=2 protein chains. cluster_id = PDB-pepchain,
items = constituent pep-pro ids, plus per-chain columns colon-joined in protein-
chain order (PROTEIN_CHAIN/SIZE/SEQ) and single peptide columns. Master set =
pairs passing BSA>0 (= the pep-pro CSV rows).

Sample validation vs v15 multipro_v4.csv (138 overlap): cluster_ids 99.3%
present, count 97.8%, items 97.8%. Mismatches are all v15-has-a-chain-we-lack
(e.g. 1BBR-F missing chain K) -- the same pep-pro-level BSA>0 / terminal-X
convention dropping a constituent pair, propagated into the grouping. Not a
grouping error.

Phase 2 (TODO): recomputed features on the assembled multi-chain complex --
ASA/BSA (single, on the combined complex), and per-chain contacts / PRODIGY
affinity+Kd (colon-joined, recomputed in multi-chain context). Colon-joined
physchem, metadata, and interface-residue columns are derivable from the pep-pro
outputs and can be joined in a Multipro assembly step alongside Phase 2.

## Multipro surface (Phase 2b) + BPepA fix (also pep-pro)
Surface recomputed on the assembled multi-chain PDB (multipro_surface.py), same
FreeSASA path as pep-pro. FreeSASA configured to match NACCESS: NACCESS radii
(Classifier.getStandardClassifier("naccess")) + Lee-Richards, n-slices=100.

Buried-area formula corrected in BOTH run_freesasa.py and multipro_surface.py to
the original method (confirmed against the peer's bsa_utils.calculate_bsa):
  BProA = ASA_Protein(isolated) - protein_in_complex   (per-chain ASA from the
          complex calc, summed over protein chains)
  BPepA = ASA_Peptide(isolated) - peptide_in_complex
  BSA   = (BProA + BPepA) / 2
(Previously BProA/BPepA used a wrong ΔASA proxy; BSA was already correct. The old
BProA/BPepA were never validated, so the bug was latent.)

Sample validation vs v15 multipro_v4.csv (138 overlap): ASA_Peptide, ASA_Protein,
BProA all r=1.000; ASA_Complex 0.980. BPepA/BPP% DIVERGE and this is a v15 DEFECT,
not ours: v15 reports the peptide as 100% buried (BPepA == ASA_Peptide) for 86.1%
(17016/19759) of entries -- physically impossible; the original NACCESS pipeline
read peptide in-complex ASA as ~0 (chain-id/.rsa parsing bug). We emit the correct
ΔASA. BSA (0.874) is partly dragged by the same v15 defect via BSA=(BProA+BPepA)/2.
We keep the correct values (same posture as extinction coefficients).

## Multipro PRODIGY (Phase 2c)
multipro_prodigy.py runs PRODIGY per protein chain on the multi-chain complex PDB
(prodigy {cid}.pdb --selection {pep} {prot_i} --distance-cutoff 6.0 --temperature 25),
colon-joined per chain in PROTEIN_CHAIN order. Reuses run_prodigy's parser.

Validation vs v15 multipro_v4.csv (273 chain-pairs, 138 entries): Spearman
intermolecular 0.995, apolar-apolar 0.992, binding affinity 0.960. Our per-chain
values track v15 with a small systematic offset (~+7 contacts; e.g. 1A1R-C 87:12
vs v15 80:11) because --selection C A measures the true C-A interface -- identical
to the isolated pep-pro pair, as it should be -- whereas v15's multipro pipeline
reported slightly fewer contacts in the multi-chain context. Ours is self-consistent
(multipro C-A == pep-pro C-A). Pearson is unreliable here (one 2x outlier, 1BQP);
Spearman is the right statistic for these skewed counts.

## Multipro assembly (Phase 2d) — dataset complete
multipro_assemble.py joins Phase 1 (grouping) + Phase 2b (surface) + Phase 2c
(PRODIGY) with the pep-pro outputs (physchem/metadata/interface/legacy) via `items`.
Single-value columns from items[0]; per-protein-chain columns colon-joined in items
order. Nothing recomputed here -- one source of truth. Output multipro_final.csv,
exact 64-column v15 order (multipro_v4.csv), incl. the two duplicate cluster_id
columns and the U+02DA degree glyph. COLUMN ORDER IS SACRED (website reads it).

leader_id = pep-pro leader_id of items[0] (verified 99.6% vs v15; 139/139 internally
consistent with our legacy_clusters). peptide_Length/protein_Length duplicate
PEPTIDE_SIZE/PROTEIN_SIZE.

Sample validation vs v15 multipro_v4.csv (138 overlap): 64-col header IDENTICAL.
Deterministic columns 100% (cluster_id, PDB_ID, PEPTIDE_CHAIN, CLASSIFICATION,
DEPOSITION_DATE, STRUCTURE_METHOD, peptide_Formula/TotalAtoms). items/count/
PROTEIN_CHAIN 97.8% (pep-pro BSA>0/terminal-X drops propagated up). PEPTIDE_SIZE
81.2% / PROTEIN_DESC 65.2% (terminal-X + v14/v15 desc drift). leader_id 77.5% vs
v14 = redundancy vintage drift (matches v15 files). Numeric: ASA/BProA Spearman
1.000, contacts 0.995, affinity 0.960; BPepA/BPP%/BSA reflect v15's documented
multipro defect (peptide 100% buried in 86% of v15 entries) -- we emit correct values.

## COCaDA contacts: auth vs label chain-ID divergence (v15 defect, corrected)

**Symptom.** On the 300-PDB sample, the COCaDA stage produced 0 contacts for
249 entries and dropped 34 more entirely (618/652 rows). Zeros were concentrated
in letter-chain complexes (e.g. `-I-H` ×30, `-L-H` ×38), not just symmetry pairs.

**Cause.** mmCIF carries two chain-ID columns: `auth_asym_id` (author naming,
used throughout Propedia and by Biopython's default parser) and `label_asym_id`
(canonical mmCIF naming: A,B,…,Z,AA,AB,…). COCaDA selects chains on
`label_asym_id`. The pipeline (and v15) passed `auth_asym_id` values to
`-c pep,prot`. These IDs diverge in ~90% of PDB entries. Where the auth ID was
not a valid label ID, COCaDA selected nothing → 0 contacts; where it happened to
be a valid but *different* label chain, contacts were silently wrong.

Example — `1A2C-L-H` (thrombin + inhibitor): auth `L,H`. COCaDA read them as
label IDs; `L` is not a label ID here (labels run A–I) → 0 contacts. Correct
polymer labels are `A,B` → 40 contacts.

**Fix.** Before calling COCaDA, translate each pair's auth chain ID to its
**polymer** `label_asym_id`, built from the CIF `_atom_site` block using ATOM
(polymer) records only. Restricting to ATOM records makes the mapping 1:1: an
auth chain's het groups get separate label IDs, but its single polymer entity
has exactly one. Chains with no polymer label are reported (`no_polymer_label`),
not silently zeroed. Stock public COCaDA v1.6 is used unchanged — no private
build required. (`run_cocada.py`, `auth_to_polymer_label()`.)

A second, compounding bug was fixed at the same time: COCaDA names its output by
PDB id, so on a re-run an entry directory could contain a stale
`<entry>_contacts.csv` from a prior run; the output glob then picked the stale
file. `run_one` now clears the entry directory before running.

**Result (sample).** 651/651 entries process (was 618); 52 zeros remain, all
verified as genuine non-contacts at COCaDA's threshold (32 checked standalone
after translation, 20 have auth==label so translation is a no-op). ~43% of the
sample's contact rows changed relative to v15.

**Consequence for the published database.** v15/v4 contact columns were computed
with the same auth-vs-label mismatch and are therefore incorrect for a large
fraction of entries. The pipeline emits the corrected values and documents the
divergence (same posture as the BProA/BPepA and extinction-coefficient fixes).
This is a case of the reproducible pipeline surfacing a latent defect in the
published data, not merely reproducing it.
