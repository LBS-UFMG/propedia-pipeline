rule ml_classifiers:
    input:
        signatures = f"{OUT}/seq_signatures.tsv",
    output:
        report = f"{OUT}/ml_report.tsv",
        best   = f"{OUT}/best_model.tsv",
    params:
        train_dir    = config["ml"]["train_dir"],
        test_dir     = config["ml"]["test_dir"],
        ifeature_dir = config["machine"]["ifeature_dir"],
        classes      = config["ml"]["classes"],
    script:
        "../scripts/run_ml.py"
		
rule therapeutic:
    input:
        pairs = f"{OUT}/pairs.tsv",
        best  = f"{OUT}/best_model.tsv",
    output:
        scores = f"{OUT}/therapeutic.tsv",
    params:
        train_dir    = config["ml"]["train_dir"],
        ifeature_dir = config["machine"]["ifeature_dir"],
    script:
        "../scripts/therapeutic_scoring.py"