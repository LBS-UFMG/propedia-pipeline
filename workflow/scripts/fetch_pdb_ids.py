"""Rule 1: candidate PDB IDs via a permissive, structure-agnostic prefilter.

An entry qualifies if, as of the snapshot's initial_release_date, it BOTH:
  (a) contains a Protein polymer entity with SEQRES length >= protein_len_min, AND
  (b) has >= min_polymer_chains deposited polymer chains (a possible partner).

We intentionally do NOT prefilter on peptide length or total monomer count here:
peptides are often short *modeled* fragments of long SEQRES chains, and the size
cap concerns amino acids (not rRNA monomers). Both are decided in Rule 3 from the
actual structure. This keeps the candidate net a safe superset of the database.
"""
import json
import sys
import time
import urllib.parse
import urllib.request

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
PAGE_ROWS = 10000


def date_node(snapshot_date):
    return {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_accession_info.initial_release_date",
        "operator": "less_or_equal", "value": snapshot_date}}


def _get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=180) as resp:
                if resp.status == 204:
                    return None
                return json.load(resp)
        except Exception as exc:                       # noqa: BLE001
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))
            print(f"  retry {attempt + 1}: {exc}", file=sys.stderr)


def search_entries(nodes):
    ids, start = set(), 0
    while True:
        request = {
            "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": start, "rows": PAGE_ROWS},
                "results_verbosity": "compact",
            },
        }
        url = f"{SEARCH_URL}?json={urllib.parse.quote(json.dumps(request))}"
        payload = _get(url)
        if not payload:
            break
        batch = payload.get("result_set", [])
        if not batch:
            break
        ids.update(x.split("_")[0] for x in batch)
        start += PAGE_ROWS
        if start >= payload.get("total_count", 0):
            break
        time.sleep(1)
    return ids


def load_ids(path):
    with open(path) as fh:
        return {line.strip().upper() for line in fh if line.strip()}


def main():
    p = snakemake.params                                   # noqa: F821
    date = [date_node(p.snapshot_date)]

    protein_entries = search_entries(date + [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "entity_poly.rcsb_entity_polymer_type",
            "operator": "exact_match", "value": p.protein_polymer_type}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "entity_poly.rcsb_sample_sequence_length",
            "operator": "greater_or_equal", "value": p.protein_len_min}},
    ])
    multichain_entries = search_entries(date + [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.deposited_polymer_entity_instance_count",
            "operator": "greater_or_equal", "value": p.min_polymer_chains}},
    ])

    candidates = sorted(protein_entries & multichain_entries)

    meta = {
        "strategy": "prefilter-permissive",
        "snapshot_date": p.snapshot_date,
        "n_protein_entries": len(protein_entries),
        "n_multichain_entries": len(multichain_entries),
        "n_candidates": len(candidates),
    }

    missing = []
    oracle_path = getattr(p, "oracle_ids", "") or ""
    if oracle_path:
        oracle = load_ids(oracle_path)
        missing = sorted(oracle - {c.upper() for c in candidates})
        meta["oracle_size"] = len(oracle)
        meta["oracle_missing"] = len(missing)
        meta["oracle_recall"] = round(1 - len(missing) / len(oracle), 4) if oracle else None

    with open(snakemake.output.ids, "w") as fh:            # noqa: F821
        fh.write("\n".join(candidates) + "\n")
    with open(snakemake.output.missing, "w") as fh:        # noqa: F821
        fh.write("\n".join(missing) + ("\n" if missing else ""))
    with open(snakemake.output.meta, "w") as fh:           # noqa: F821
        json.dump(meta, fh, indent=2)

    print(f"candidates: {len(candidates)}  oracle_missing: {len(missing)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()