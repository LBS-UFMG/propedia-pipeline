"""Rule 3 (sample): from each candidate CIF, classify chains by MODELED
amino-acid count, find peptide(2-50)-protein(>50) pairs in contact (<=cutoff A),
write the two-chain complex as a clean ATOM-only mmCIF, and record sequences.

mmCIF (not PDB) is the intermediate: it is lossless for large complexes
(>99,999 atoms), multi-character author chain IDs, and wide residue numbering,
none of which legacy PDB can represent. Biopython MMCIFIO preserves the author
chain IDs (auth_asym_id == chain.id); every downstream consumer keys on those
same IDs (the pep_chain/prot_chain recorded here), and both PRODIGY and FreeSASA
read auth_asym_id from the written file, so chain selection stays unambiguous.
"""
import gzip
import os
import sys
import tempfile
from types import SimpleNamespace

from Bio.PDB import MMCIFParser, MMCIFIO, NeighborSearch, Select

import checkpoint

VERSION = "3"            # bump on logic change -> invalidates old checkpoints
                         # (v3: v15-style X padding — polymer_set residue membership)
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1_extended


def one(resname):
    code = protein_letters_3to1_extended.get(resname.capitalize(), "X")
    return code if len(code) == 1 else "X"


def polymer_set(cif_text):
    """(auth_asym_id, auth_seq_id) of every POLYMER residue = `_atom_site` rows whose
    `label_seq_id` is numeric. These are the peptide/protein residues INCLUDING terminal
    caps (ACE/NH2) and modified residues; waters/ligands/ions carry '.'/'?' and are
    excluded. This is what restores v15's convention: keep every polymer residue in the
    sequence, writing non-standard ones as `X` (see `one`). Returns an EMPTY set if the
    cif has no `label_seq_id` column, which makes `modeled_aa` fall back to the
    amino-acid test (older/degenerate files never regress)."""
    poly, cols, started = set(), {}, False
    for line in cif_text.splitlines():
        s = line.strip()
        if s.startswith("_atom_site."):
            cols[s.split(".", 1)[1]] = len(cols)
        elif cols and (s.startswith("ATOM ") or s.startswith("HETATM")):
            started = True
            f = s.split()
            try:
                if f[cols["label_seq_id"]] not in (".", "?"):
                    poly.add((f[cols["auth_asym_id"]], f[cols["auth_seq_id"]]))
            except (KeyError, IndexError):
                return set()                     # unexpected layout -> fall back
        elif started:
            break                                # past the _atom_site loop
    return poly


def modeled_aa(chain, poly=None):
    """Residues that make up the peptide/protein sequence. v15 convention: every
    POLYMER residue (standard, modified, or capped) counts — non-standard ones become
    `X` in `seq_of`. `poly` is the polymer-residue set from `polymer_set`; when it is
    empty/None (no label_seq_id in the cif) we fall back to the amino-acid test."""
    if poly:
        cid = chain.id
        return [r for r in chain if (cid, str(r.id[1])) in poly]
    return [r for r in chain if is_aa(r, standard=False)]


def seq_of(residues):
    return "".join(one(r.resname) for r in residues)


class PairSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)

    def accept_model(self, model):
        return True

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

    def accept_residue(self, residue):
        return is_aa(residue, standard=False)


def shard_path(cif_dir, pid):
    sub = pid[1:3].lower() if len(pid) >= 3 else "_"
    return os.path.join(cif_dir, sub, f"{pid}.cif.gz")


def count_atoms(cif_dir, pid):
    """Cheap atom count: scan the gzipped CIF text for _atom_site ATOM/HETATM rows,
    WITHOUT a Biopython parse (which is the expensive step we want to avoid for huge
    structures). Returns an int, or None if the file is missing/unreadable."""
    path = shard_path(cif_dir, pid)
    try:
        n = 0
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("ATOM ") or line.startswith("HETATM"):
                    n += 1
        return n
    except Exception:                                          # noqa: BLE001
        return None


def load_atom_cache(path):
    """pid -> atom count, from the durable cache (empty dict if absent)."""
    cache = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            next(fh, None)                                      # header
            for line in fh:
                pid, _, n = line.rstrip("\n").partition("\t")
                if pid and n.isdigit():
                    cache[pid] = int(n)
    return cache


def save_atom_cache(path, cache):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write("pid\tn_atoms\n")
        for pid in sorted(cache):
            fh.write(f"{pid}\t{cache[pid]}\n")
    os.replace(tmp, path)


def load_first_model(pid, cif_dir, parser):
    path = shard_path(cif_dir, pid)
    with gzip.open(path, "rt") as fh:
        data = fh.read()
    poly = polymer_set(data)                # v15-style residue membership (label_seq_id)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False)
    try:
        tmp.write(data)
        tmp.close()
        structure = parser.get_structure(pid, tmp.name)
    finally:
        os.unlink(tmp.name)
    return next(iter(structure)), poly      # first model (first NMR conformer)


def process(pid, p, io, parser):
    rows = []
    try:
        model, poly = load_first_model(pid, p.cif_dir, parser)
    except Exception as exc:                # noqa: BLE001
        return rows, f"parse_error"
    chain_res = {c.id: modeled_aa(c, poly) for c in model}
    chain_res = {c: r for c, r in chain_res.items() if r}
    peptides = {c: r for c, r in chain_res.items()
                if p.pep_min <= len(r) <= p.pep_max}
    proteins = {c: r for c, r in chain_res.items() if len(r) >= p.prot_min}
    if not peptides or not proteins:
        return rows, "no_pair"

    ns = NeighborSearch(list(model.get_atoms()))
    for pep_id, pep_res in peptides.items():
        for prot_id, prot_res in proteins.items():
            if pep_id == prot_id:
                continue
            contact = False
            for atom in (a for r in pep_res for a in r.get_atoms()):
                for near in ns.search(atom.coord, p.cutoff, level="A"):
                    if near.get_parent().get_parent().id == prot_id:
                        contact = True
                        break
                if contact:
                    break
            if not contact:
                continue
            entry_id = f"{pid}-{pep_id}-{prot_id}"
            io.set_structure(model)
            io.save(os.path.join(p.cif_out_dir, f"{entry_id}.cif"),
                    PairSelect([prot_id, pep_id]))
            rows.append({"id": entry_id, "pdb": pid,
                         "pep_chain": pep_id, "prot_chain": prot_id,
                         "pep_size": len(pep_res), "prot_size": len(prot_res),
                         "pep_seq": seq_of(pep_res),
                         "prot_seq": seq_of(prot_res)})
    return rows, "ok"


def extract_worker(item):
    """worker(item) -> (rows|[], status). Parses one input CIF, writes that entry's
    per-pair CIFs (side effect), and returns the pair rows. Deterministic, so every
    outcome is checkpointed (never retried)."""
    pid, cif_dir, cif_out_dir, pep_min, pep_max, prot_min, cutoff = item
    p = SimpleNamespace(cif_dir=cif_dir, cif_out_dir=cif_out_dir, pep_min=pep_min,
                        pep_max=pep_max, prot_min=prot_min, cutoff=cutoff)
    rows, status = process(pid, p, MMCIFIO(), MMCIFParser(QUIET=True))
    return rows, status


def main():
    p = snakemake.params                                   # noqa: F821
    threads = getattr(snakemake, "threads", 1) or 1        # noqa: F821
    os.makedirs(p.cif_out_dir, exist_ok=True)              # resume: do NOT clear
    sample = [l.strip() for l in open(snakemake.input.ids) if l.strip()]  # noqa: F821

    # ---- size gate (membership filter, NOT a compute/checkpoint param) ----
    # Skip structures with more than `max_atoms` atoms so huge complexes (ribosomes,
    # etc.) never enter pairs.tsv and never cascade downstream. Atom counts are cached
    # durably so the decision is instant on re-runs and only NEW pids are scanned;
    # raising the limit later simply re-admits previously-skipped pids as new work.
    max_atoms = int(getattr(p, "max_atoms", 0) or 0)
    atom_cache_path = getattr(p, "atom_cache", None)
    cache = load_atom_cache(atom_cache_path) if atom_cache_path else {}
    scanned = 0
    if max_atoms > 0:
        for pid in sample:
            if pid not in cache:
                n = count_atoms(p.cif_dir, pid)
                if n is not None:
                    cache[pid] = n
                    scanned += 1
        if atom_cache_path and scanned:
            save_atom_cache(atom_cache_path, cache)
    included, oversized = [], []
    for pid in sample:
        n = cache.get(pid)
        if max_atoms > 0 and n is not None and n > max_atoms:
            oversized.append((pid, n))
        else:
            included.append(pid)
    if max_atoms > 0:
        print(f"size gate: max_atoms={max_atoms}; {len(included)} kept, "
              f"{len(oversized)} skipped ({scanned} newly scanned)", file=sys.stderr)
        over_path = getattr(p, "oversized", None)
        if over_path:
            os.makedirs(os.path.dirname(over_path) or ".", exist_ok=True)
            with open(over_path, "w") as fh:
                fh.write("pdb\tn_atoms\n")
                for pid, n in sorted(oversized):
                    fh.write(f"{pid}\t{n}\n")

    items = [(pid, p.cif_dir, p.cif_out_dir, p.pep_min, p.pep_max, p.prot_min,
              p.cutoff) for pid in included]
    workdir = checkpoint.namespace(p.ckpt, VERSION,
                                   {"pep_min": p.pep_min, "pep_max": p.pep_max,
                                    "prot_min": p.prot_min, "cutoff": p.cutoff})
    results = checkpoint.run(items, extract_worker, workdir, threads=threads,
                             id_of=lambda it: it[0], stage="extract")

    cols = ["id", "pdb", "pep_chain", "prot_chain",
            "pep_size", "prot_size", "pep_seq", "prot_seq"]
    all_rows, keep_ids = [], set()
    for pid in included:
        res = results.get(pid)
        if not res or not res["record"]:
            continue
        for r in res["record"]:
            all_rows.append(r)
            keep_ids.add(r["id"])
    with open(snakemake.output.pairs, "w") as fh:          # noqa: F821
        fh.write("\t".join(cols) + "\n")
        for r in all_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # hygiene: prune per-pair CIFs no longer referenced (smaller/different sample);
    # replaces the old start-of-run clear, which would have wiped resume progress.
    for f in os.listdir(p.cif_out_dir):
        if f.endswith(".cif") and f[:-4] not in keep_ids:
            os.remove(os.path.join(p.cif_out_dir, f))

    print(f"DONE pairs={len(all_rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()