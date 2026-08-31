"""Multipro Phase 2a: write the multi-chain mmCIF for each Multipro entry
(peptide chain + all its protein chains), from the source CIF. Reuses extract_pairs
for CIF parsing / chain selection so atom content matches the pep-pro files. The
model is cached per PDB (several Multipro entries can share one structure).

mmCIF output (not PDB): a peptide bound by several protein chains can exceed PDB's
atom/chain-ID limits; MMCIFIO preserves author chain IDs, which the surface and
PRODIGY stages select on.
"""
import csv
import os
import sys

from Bio.PDB import MMCIFParser, MMCIFIO

import extract_pairs as ep


def main():
    p = snakemake.params                                   # noqa: F821
    os.makedirs(p.cif_out_dir, exist_ok=True)              # resume: do NOT clear
    rows = list(csv.DictReader(open(snakemake.input.multipro), delimiter="\t"))  # noqa: F821
    parser = MMCIFParser(QUIET=True)
    io = MMCIFIO()
    cache = {}
    n = skipped = 0
    want = {".written"}
    for r in rows:
        cid = r["cluster_id"]
        out = os.path.join(p.cif_out_dir, f"{cid}.cif")
        want.add(f"{cid}.cif")
        if os.path.exists(out):        # resume: already written on a prior run
            skipped += 1
            continue
        pid = r["PDB_ID"]
        chains = {r["PEPTIDE_CHAIN"]} | set(r["PROTEIN_CHAIN"].split(":"))
        if pid not in cache:
            try:
                cache[pid] = ep.load_first_model(pid, p.cif_dir, parser)[0]  # model only
            except Exception:                              # noqa: BLE001
                cache[pid] = None
        model = cache[pid]
        if model is None:
            want.discard(f"{cid}.cif")
            continue
        io.set_structure(model)
        tmp = out + ".tmp"
        io.save(tmp, ep.PairSelect(chains))
        os.replace(tmp, out)           # atomic: no partial CIF on a kill
        n += 1
        if n % 50 == 0:
            print(f"{n} multipro cifs ({skipped} resumed)", file=sys.stderr)

    for f in os.listdir(p.cif_out_dir):     # prune stale (smaller/different sample)
        if f.endswith(".cif") and f not in want:
            os.remove(os.path.join(p.cif_out_dir, f))

    open(os.path.join(p.cif_out_dir, ".written"), "w").write(f"{n + skipped}\n")
    print(f"DONE {n} new, {skipped} resumed multipro cifs", file=sys.stderr)


if __name__ == "__main__":
    main()
