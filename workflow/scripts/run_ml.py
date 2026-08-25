"""Bake-off for the six therapeutic-peptide classifiers (AAP/ABP/ACP/AIP/QSP/SBP).
All six models are trained per class and scored on the held-out test set; the
per-class winner (by test-set BRIER SCORE, i.e. probability quality) is exported
for the scoring stage. Brier, not AUC, drives selection because the shipped
columns are probabilities: a model can rank well (high AUC) yet emit poorly
calibrated probabilities (e.g. SVM's Platt scaling). AUC is still reported.

Data: labeled peptide SEQUENCES from the propedia26-sm train/test TSVs (col 1 =
peptide, col 2 = target). Features are RECOMPUTED via the same iFeature 9-descriptor
1248-vector used by the signature and scoring stages, so the bake-off ranking and
the shipped scores live in one feature space. StandardScaler sits inside each
pipeline (fit on train only -> no leakage; a no-op for trees, helps SVM/NN).
Default hyperparameters (proof-of-principle, not optimized)."""
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
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (accuracy_score, roc_auc_score, f1_score,
                             matthews_corrcoef, brier_score_loss)

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


def read_tsv(path):
    """(clean_seq, label01) from a propedia26-sm TSV (col1 seq, col2 target)."""
    out = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            seq = "".join(c for c in f[0].upper() if c in STD_AA)
            if len(seq) < 3:
                continue
            out.append((seq, 1 if f[1].strip().lower().startswith("pos") else 0))
    return out


def ifeature_matrix(seqs_by_id, ifeature_dir):
    tmp = tempfile.mkdtemp()
    fasta = os.path.join(tmp, "seqs.fasta")
    with open(fasta, "w") as fh:
        for sid, seq in seqs_by_id.items():
            fh.write(f">{sid}\n{seq}\n")
    feats = {}
    for desc in CORE_9:
        out_tsv = os.path.join(tmp, f"{desc}.tsv")
        subprocess.run(["python3", "iFeature.py", "--file", os.path.abspath(fasta),
                        "--type", desc, "--out", os.path.abspath(out_tsv)],
                       cwd=ifeature_dir, capture_output=True, text=True)
        with open(out_tsv) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                feats.setdefault(f[0], []).extend(f[1:])
    return {sid: np.array(v, dtype=float) for sid, v in feats.items() if len(v) == 1248}


def run_class(cls, train_dir, test_dir, ifeature_dir):
    train = read_tsv(os.path.join(train_dir, f"train_{cls}.tsv"))
    test = read_tsv(os.path.join(test_dir, f"test_main_{cls}.tsv"))
    tr = {f"tr_{i}": s for i, (s, _) in enumerate(train)}
    te = {f"te_{i}": s for i, (s, _) in enumerate(test)}
    allf = ifeature_matrix({**tr, **te}, ifeature_dir)
    Xtr = np.array([allf[f"tr_{i}"] for i, _ in enumerate(train) if f"tr_{i}" in allf])
    ytr = np.array([y for i, (_, y) in enumerate(train) if f"tr_{i}" in allf])
    Xte = np.array([allf[f"te_{i}"] for i, _ in enumerate(test) if f"te_{i}" in allf])
    yte = np.array([y for i, (_, y) in enumerate(test) if f"te_{i}" in allf])

    results = []
    for name, make in MODELS.items():
        try:
            m = make_pipeline(StandardScaler(), make()).fit(Xtr, ytr)
            pred = m.predict(Xte)
            proba = m.predict_proba(Xte)[:, 1]
            auc = roc_auc_score(yte, proba) if len(set(yte)) > 1 else float("nan")
            brier = brier_score_loss(yte, proba)
            results.append((name, accuracy_score(yte, pred), auc,
                            f1_score(yte, pred), matthews_corrcoef(yte, pred), brier))
        except Exception as exc:                       # noqa: BLE001
            results.append((name, float("nan"), float("nan"), float("nan"),
                            float("nan"), float("nan")))
            print(f"  {name}: ERROR {exc}", file=sys.stderr)
    return len(Xtr), len(Xte), results


def main():
    p = snakemake.params                                   # noqa: F821
    train_dir = os.path.expanduser(p.train_dir)
    test_dir = os.path.expanduser(p.test_dir)
    ifeature_dir = os.path.expanduser(p.ifeature_dir)
    classes = p.classes

    best = {}
    with open(snakemake.output.report, "w") as out:        # noqa: F821
        out.write("class\tn_train\tn_test\tmodel\taccuracy\tauc\tf1\tmcc\tbrier\n")
        for cls in classes:
            print(f"=== {cls} ===", file=sys.stderr)
            ntr, nte, res = run_class(cls, train_dir, test_dir, ifeature_dir)
            best_brier, best_model, best_auc = float("inf"), None, float("nan")
            for name, acc, auc, f1, mcc, brier in res:
                out.write(f"{cls}\t{ntr}\t{nte}\t{name}\t{acc:.3f}\t{auc:.3f}\t"
                          f"{f1:.3f}\t{mcc:.3f}\t{brier:.4f}\n")
                print(f"  {name:20s} auc={auc:.3f} brier={brier:.4f}", file=sys.stderr)
                if brier == brier and brier < best_brier:   # brier==brier filters NaN
                    best_brier, best_model, best_auc = brier, name, auc
            best[cls] = (best_model, best_brier, best_auc)

    with open(snakemake.output.best, "w") as bf:           # noqa: F821
        bf.write("class\tmodel\tbrier\tauc\n")
        for cls in classes:
            m, br, a = best[cls]
            bf.write(f"{cls}\t{m}\t{br:.4f}\t{a:.3f}\n")
    print("DONE", file=sys.stderr)


if __name__ == "__main__":
    main()