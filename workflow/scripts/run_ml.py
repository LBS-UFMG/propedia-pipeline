"""Reproduce the six therapeutic-peptide classifiers (AAP/ABP/ACP/AIP/QSP/SBP).
Features: iFeature 1248-vector (9 descriptors, matching tonight's signature stage).
Models: the Orange bake-off ported to scikit-learn (SVM, GradientBoosting,
LogisticRegression, kNN, NaiveBayes, MLP), scored on the held-out test set.
Reproduces the METHOD; matching exact published numbers needs Suppl. Tables S1-S7."""
import csv
import os
import subprocess
import sys
import tempfile

import numpy as np
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, matthews_corrcoef

CORE_9 = ["AAC", "DPC", "DDE", "GAAC", "GDPC", "GTPC", "CTDC", "CTDT", "CTDD"]
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")

MODELS = {
    "SVM": lambda: SVC(probability=True, random_state=0),
    "GradientBoosting": lambda: GradientBoostingClassifier(random_state=0),
    "LogisticRegression": lambda: LogisticRegression(max_iter=2000, random_state=0),
    "kNN": lambda: KNeighborsClassifier(),
    "NaiveBayes": lambda: GaussianNB(),
    "NeuralNet": lambda: MLPClassifier(max_iter=1000, random_state=0),
}


def read_pairs(path):
    """Return list of (id, peptide, label01) from a cs4 5-column file."""
    out = []
    for r in csv.DictReader(open(path)):
        seq = "".join(c for c in r["Peptide"].upper() if c in STD_AA)
        if len(seq) < 3:
            continue
        label = 1 if r["Actual"].strip().lower().startswith("pos") else 0
        out.append((r["ID"], seq, label))
    return out


def ifeature_matrix(seqs_by_id, ifeature_dir):
    """Run iFeature over a FASTA of sequences; return {id: [floats]} for the 1248 vec."""
    tmp = tempfile.mkdtemp()
    fasta = os.path.join(tmp, "seqs.fasta")
    with open(fasta, "w") as fh:
        for sid, seq in seqs_by_id.items():
            fh.write(f">{sid}\n{seq}\n")
    feats = {}
    order = None
    for desc in CORE_9:
        out_tsv = os.path.join(tmp, f"{desc}.tsv")
        subprocess.run(["python3", "iFeature.py", "--file", os.path.abspath(fasta),
                        "--type", desc, "--out", os.path.abspath(out_tsv)],
                       cwd=ifeature_dir, capture_output=True, text=True)
        with open(out_tsv) as fh:
            fh.readline()  # header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                feats.setdefault(f[0], []).extend(f[1:])
    return {sid: np.array(v, dtype=float) for sid, v in feats.items()
            if len(v) == 1248}


def run_class(cls, csv_dir, ifeature_dir):
    train = read_pairs(os.path.join(csv_dir, f"{cls}_train.csv"))
    test = read_pairs(os.path.join(csv_dir, f"{cls}_test_main.csv"))
    # unique ids per split (prefix to avoid train/test id collision)
    tr_seqs = {f"tr_{i}": s for i, (_, s, _) in enumerate(train)}
    te_seqs = {f"te_{i}": s for i, (_, s, _) in enumerate(test)}
    allf = ifeature_matrix({**tr_seqs, **te_seqs}, ifeature_dir)

    Xtr, ytr, Xte, yte = [], [], [], []
    for i, (_, s, y) in enumerate(train):
        k = f"tr_{i}"
        if k in allf:
            Xtr.append(allf[k]); ytr.append(y)
    for i, (_, s, y) in enumerate(test):
        k = f"te_{i}"
        if k in allf:
            Xte.append(allf[k]); yte.append(y)
    Xtr, ytr, Xte, yte = map(np.array, (Xtr, ytr, Xte, yte))

    results = []
    for name, make in MODELS.items():
        try:
            m = make().fit(Xtr, ytr)
            pred = m.predict(Xte)
            proba = m.predict_proba(Xte)[:, 1] if hasattr(m, "predict_proba") else pred
            results.append((name,
                            accuracy_score(yte, pred),
                            roc_auc_score(yte, proba) if len(set(yte)) > 1 else float("nan"),
                            f1_score(yte, pred),
                            matthews_corrcoef(yte, pred)))
        except Exception as exc:                       # noqa: BLE001
            results.append((name, float("nan"), float("nan"), float("nan"), str(exc)[:30]))
    return len(Xtr), len(Xte), results


def main():
    p = snakemake.params                                   # noqa: F821
    csv_dir = os.path.expanduser(p.csv_dir)
    ifeature_dir = os.path.expanduser(p.ifeature_dir)
    classes = p.classes

    with open(snakemake.output.report, "w") as out:        # noqa: F821
        out.write("class\tn_train\tn_test\tmodel\taccuracy\tauc\tf1\tmcc\n")
        for cls in classes:
            print(f"=== {cls} ===", file=sys.stderr)
            ntr, nte, res = run_class(cls, csv_dir, ifeature_dir)
            for name, acc, auc, f1, mcc in res:
                out.write(f"{cls}\t{ntr}\t{nte}\t{name}\t{acc:.3f}\t"
                          f"{auc if isinstance(auc,str) else f'{auc:.3f}'}\t"
                          f"{f1 if isinstance(f1,str) else f'{f1:.3f}'}\t"
                          f"{mcc if isinstance(mcc,str) else f'{mcc:.3f}'}\n")
                print(f"  {name:20s} acc={acc:.3f} auc={auc:.3f} f1={f1:.3f}",
                      file=sys.stderr)
    print("DONE", file=sys.stderr)


if __name__ == "__main__":
    main()