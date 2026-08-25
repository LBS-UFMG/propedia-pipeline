"""Therapeutic-class scoring: six per-peptide probabilities (AAP, ABP, ACP, AIP,
QSP, SBP) from GradientBoosting classifiers.

Training peptides (labeled sequences from train_{CLASS}.tsv) and Propedia peptides
are featurized through the SAME iFeature pass (reusing run_ifeature) so train and
predict share one feature space. Per class: StandardScaler + GradientBoosting with
default hyperparameters (proof-of-principle, not optimized). Output: id + 6
probabilities (2 decimals), keyed by peptide so all pairs sharing a peptide match.
"""
import csv
import os
import sys
import tempfile

import numpy as np
import run_ml as ml_models   # reuse the MODELS registry (one source of truth)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import run_ifeature as ifx   # reuse CORE_9 / STD_AA / run_descriptor / load_matrix

CLASSES = ["AAP", "ABP", "ACP", "AIP", "QSP", "SBP"]
RANDOM_STATE = 42


def clean(seq):
    return "".join(c for c in seq.upper() if c in ifx.STD_AA)


def load_train_labels(train_dir):
    """{class: {clean_seq: 1/0}} from train_{CLASS}.tsv (col1 seq, col2 target)."""
    out = {}
    for cls in CLASSES:
        path = os.path.join(train_dir, f"train_{cls}.tsv")
        d = {}
        with open(path) as fh:
            next(fh)  # header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                seq = clean(f[0])
                if len(seq) < 3:
                    continue
                d[seq] = 1 if f[1].strip().lower().startswith("pos") else 0
        out[cls] = d
    return out


def featurize(seqs, ifeature_dir):
    """{seq: np.array} via the same 9 iFeature descriptors as run_ifeature."""
    tmp = tempfile.mkdtemp()
    fasta = os.path.join(tmp, "seqs.fasta")
    with open(fasta, "w") as fh:
        for s in seqs:
            fh.write(f">{s}\n{s}\n")     # header == sequence (unique key)
    headers, feats, order = [], {}, None
    for desc in ifx.CORE_9:
        out_tsv = os.path.join(tmp, f"{desc}.tsv")
        if not ifx.run_descriptor(ifeature_dir, os.path.abspath(fasta),
                                  desc, os.path.abspath(out_tsv)):
            print(f"  {desc}: FAILED", file=sys.stderr)
            continue
        hdr, rows = ifx.load_matrix(out_tsv)
        headers += hdr
        if order is None:
            order = list(rows.keys())
        for sid, vals in rows.items():
            feats.setdefault(sid, []).extend(vals)
    width = len(headers)
    return {s: np.array(v, dtype=float) for s, v in feats.items() if len(v) == width}


def main():
    p = snakemake.params                                   # noqa: F821
    train_dir = os.path.expanduser(p.train_dir)
    ifeature_dir = os.path.expanduser(p.ifeature_dir)
    labels = load_train_labels(train_dir)

    best = {}
    with open(snakemake.input.best) as fh:                 # noqa: F821
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            best[parts[0]] = parts[1]

    pairs = list(csv.DictReader(open(snakemake.input.pairs), delimiter="\t"))  # noqa: F821
    pep_of = {r["id"]: clean(r["pep_seq"]) for r in pairs}
    prop_seqs = sorted({s for s in pep_of.values() if len(s) >= 3})
    train_seqs = {s for d in labels.values() for s in d}
    all_seqs = sorted(set(train_seqs) | set(prop_seqs))
    print(f"featurizing {len(all_seqs)} seqs ({len(train_seqs)} train, "
          f"{len(prop_seqs)} propedia)", file=sys.stderr)
    feat = featurize(all_seqs, ifeature_dir)

    scores = {s: {} for s in prop_seqs}
    pred_seqs = [s for s in prop_seqs if s in feat]
    X_pred = np.vstack([feat[s] for s in pred_seqs]) if pred_seqs else None
    for cls in CLASSES:
        items = [(s, y) for s, y in labels[cls].items() if s in feat]
        X = np.vstack([feat[s] for s, _ in items])
        y = np.array([y for _, y in items])
        model = make_pipeline(StandardScaler(), ml_models.MODELS[best[cls]]())
        model.fit(X, y)
        if X_pred is not None:
            for s, pr in zip(pred_seqs, model.predict_proba(X_pred)[:, 1]):
                scores[s][cls] = round(float(pr), 2)
        print(f"  {cls}: trained on {len(y)} (pos={int(y.sum())})", file=sys.stderr)

    with open(snakemake.output.scores, "w") as fh:         # noqa: F821
        fh.write("id\t" + "\t".join(CLASSES) + "\n")
        for r in pairs:
            s = pep_of[r["id"]]
            vals = [str(scores[s][cls]) if cls in scores.get(s, {}) else ""
                    for cls in CLASSES]
            fh.write(r["id"] + "\t" + "\t".join(vals) + "\n")
    print("DONE", file=sys.stderr)


if __name__ == "__main__":
    main()