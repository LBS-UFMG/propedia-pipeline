# Publishing to the propedia26 website — data layout & schema map

Reference for the future **`package`/publish step** that turns this pipeline's outputs
(`results/<mode>/` + `state/<mode>/`) into the file tree the [propedia26](https://github.com/LBS-UFMG/propedia26)
web app serves. Written from reading the site's `app/Controllers/Entry.php` / `Export.php`
and the committed example files. **Some column orders are marked TO-VERIFY** — confirm
against a live v17 file before relying on them.

> Status: **not implemented.** The pipeline produces the ingredients; nothing yet writes
> the site layout. This doc is the spec. It couples with the **PISA** work (see gaps).

## The website's data tree

Served from `propedia26/public/data/`. `<mode>` is the dataset (`propedia` for pep–pro;
multipro lives under a `multipro/` subdir). Paths taken verbatim from the controllers:

| Path | Delimiter | Header? | Purpose | Our source |
|---|---|---|---|---|
| `data/<mode>/csv/<id[0]>/<id>.csv` | `;` | no | **per-entry full row** (Entry page) | `results/<mode>/propedia.csv` (split) |
| `data/<mode>/contacts/<id>/<PDB>_contacts.csv` | `,` | yes | **per-entry COCaDA contacts** | `state/<mode>/cocada/<id>/<id>_contacts.csv` |
| `data/<mode>/multipro/csv/<id[0]>/<id>.csv` | `;` | no | per-entry multipro row | `results/<mode>/multipro_final.csv` (split) |
| `data/<mode>/multipro/contacts/<id>/<PDB>_contacts.csv` | `,` | yes | multipro per-entry contacts | `state/<mode>/cocada` (multipro run) |
| `data/propedia26_v17.tsv` | **TAB** | no | **Explore-page summary** (~25 cols; col 22 = PISA CSS) | derive from `propedia.csv` (subset) |
| `data/clusters/*.tsv` | TAB | — | cluster tables | our legacy/CNR/therapeutic cluster files |
| `data/pdb/<id[0]>/<id>/data.cif` | — | — | structure for PyMOL export | fetched on-demand from RCSB by the site (not a bulk deliverable) |

`id` = `PDB-pepChain-protChain` (e.g. `1WRZ-B-A`); `<PDB>` = `substr(id,0,4)`.

## Artifact 1 — per-entry full CSV (`csv/<id[0]>/<id>.csv`)

`;`-delimited, **no header**, one row = the entry. The live site (**v17**) row has **93
fields** in a DIFFERENT order than our v15-based `propedia.csv` (74 fields). Confirmed
prefix order from `Entry.php` (indices 0–52):

```
0  id            1 PDB_ID       2 TITLE        3 RESOLUTION   4 CLASSIFICATION
5  DEPOSITION_DATE 6 STRUCTURE_METHOD 7 PROTEIN_CHAIN 8 PEPTIDE_CHAIN 9 PROTEIN_SIZE
10 PEPTIDE_SIZE  11 PROTEIN_DESC 12 PEPTIDE_DESC 13 PROTEIN_SEQ  14 PEPTIDE_SEQ
15 leader_id     16 is_leader
17..34  peptide_{Length,MW,pI,InstabilityIndex,AliphaticIndex,GRAVY,HydrophobicPercent,
        PositiveResidues,NegativeResidues, C,H,N,O,S, Formula,TotalAtoms,
        ExtCoeff_Disulfide,ExtCoeff_NoDisulfide}
35..52  protein_{same 18 fields}
53..92  TO-VERIFY — contacts counts, PRODIGY energy, surface (ASA/BSA), therapeutic
        scores, cluster labels, and the PISA block (see gaps). Enumerate from a live v17
        per-entry csv before implementing.
```

**Transforms the packager must apply** (vs our `propedia.csv`):
- **Reorder** columns to the v17 order above (metadata → physchem → …), not v15's order.
- **Add per-element atom counts** `peptide_C/H/N/O/S` and `protein_C/H/N/O/S` — we don't
  emit these but they're trivially derivable (we already sum them in `physchem.py`; only
  `Formula`/`TotalAtoms` are currently output).
- **Add the PISA block** — not produced yet (deferred PISA work). Until then, per-entry
  csvs cannot match v17's 93-column layout and the site's fixed indices.
- Split by `id[0]` shard; write `;`-delimited, no header.

## Artifact 2 — per-entry COCaDA contacts (`contacts/<id>/<PDB>_contacts.csv`)

**Format is byte-for-byte identical to our `run_cocada` output** — verified:
```
Chain1,Res1,ResName1,Atom1,Chain2,Res2,ResName2,Atom2,Distance,Type
A,7,E,OE1,B,303,R,NE,2.71,SB
```
Same columns, comma delimiter, same `Type` codes (`SB`, `HB`, `HY`, …). **No reformatting.**
Packager only needs to:
- take `state/<mode>/cocada/<id>/<id>_contacts.csv`,
- write it to `data/<mode>/contacts/<id>/<PDB>_contacts.csv` (rename `<id>` → `<PDB>`).

This is COCaDA's real deliverable — so the multi-char COCaDA fix (recovering ~16% of large
structures) is load-bearing for the site, not cosmetic.

## Artifact 3 — Explore summary TSV (`propedia26_v17.tsv`)

TAB-delimited, **no header**, ~25 cols per entry — a compact table for the Explore/search
listing. Column **22** (TAB-separated) holds the **PISA CSS** string the Entry page reads
(`getPisaCss`). Approximate v17 layout (TO-VERIFY, from one sampled row `1A0N-A-B`):
`id, protein_size, peptide_size, peptide_seq, TITLE, CLASSIFICATION, is_leader, leader_id,
…, method, …, [16–21] six therapeutic scores, [22] PISA CSS, [23] dG, …`.
Derivable as a projection of `propedia.csv` **once PISA exists**.

## Artifact 4 — cluster tables (`data/clusters/*.tsv`)

`AAP/ABP/ACP/AIP/QSP/SBP.tsv`, `binding.tsv`, `interface.tsv`, `sequence.tsv`,
`redundant.tsv`, `seq100_clusters-NR.tsv`, `pdb_classes.tsv`. We already have the inputs
(the inherited legacy clusters in `clusters_v15/`, our CNR `seq100`, and the therapeutic
class assignments); the packager copies/derives these into the site layout. Confirm each
file's exact column format against the repo copies.

## What we already produce vs the gaps

| Site artifact | Ready? | Gap |
|---|---|---|
| per-entry COCaDA contacts | ✅ | just rename+place |
| cluster tables | ✅ (mostly) | confirm per-file format; new entries lack legacy labels (issue #2) |
| per-entry full CSV (v17, 93 col) | ⚠️ | reorder + per-element atoms + **PISA** + `;`/no-header split |
| Explore summary TSV | ⚠️ | projection + **PISA CSS** (col 22) |
| structures (`data.cif`) | n/a | site fetches from RCSB on demand |

## Open items (do before/with the packager)
1. **PISA integration** — produce the interface block (ΔG, interface area, H-bonds, salt
   bridges, buried atoms, the CSS string). This unblocks both the per-entry v17 csv tail
   and the Explore-tsv col 22. (Owner: user + colleague; "later" work.)
2. **Enumerate the exact v17 per-entry order for indices 53–92** and the 25-col Explore
   order from live files (the `examples/` copies are stale v15).
3. **Emit per-element atom counts** in `physchem.py` (cheap; already summed internally).
4. **Legacy clusters for new entries** (issue #2) so cluster tables are complete on updates.
5. Then write the `package` rule: `results/` + `state/` → `public/data/<mode>/…`, and have
   `release` (or a new target) produce it.
