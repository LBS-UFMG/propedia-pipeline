"""Provenance: record the release in which each entry FIRST appeared in Propedia.

A durable ledger (state/<mode>/provenance.tsv) maps entry id -> first_release
(the snapshot date of the run when it first entered the database). On an update,
existing entries keep their original first_release; only new entries get the
current snapshot date. The ledger is DURABLE state — preserve it between releases
(if it is lost, every current entry is (re)stamped with the current date, losing
history). Emits the current entries' first_release for assembly.

'Entries' = the pairs that make it into propedia.csv, i.e. those passing BSA>0
(present in surface.tsv), matching the assemble master set.
"""
import csv
import os
import sys


def main():
    p = snakemake.params                                   # noqa: F821
    date = str(p.snapshot_date)
    ledger_path = p.ledger

    ids = [r["id"] for r in csv.DictReader(open(snakemake.input.surface),  # noqa: F821
                                           delimiter="\t")]

    ledger = {}
    if os.path.exists(ledger_path):
        with open(ledger_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ledger[r["id"]] = r["FIRST_RELEASE"]

    added = 0
    for eid in ids:
        if eid not in ledger:
            ledger[eid] = date
            added += 1

    # write the durable all-time ledger atomically
    os.makedirs(os.path.dirname(ledger_path) or ".", exist_ok=True)
    tmp = ledger_path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("id\tFIRST_RELEASE\n")
        for eid in sorted(ledger):
            fh.write(f"{eid}\t{ledger[eid]}\n")
    os.replace(tmp, ledger_path)

    # emit first_release for the current entries (for assembly)
    with open(snakemake.output.provenance, "w") as fh:     # noqa: F821
        fh.write("id\tFIRST_RELEASE\n")
        for eid in ids:
            fh.write(f"{eid}\t{ledger[eid]}\n")

    print(f"provenance: {len(ids)} entries, {added} new this release "
          f"(first_release={date}); ledger holds {len(ledger)} all-time",
          file=sys.stderr)


if __name__ == "__main__":
    main()
