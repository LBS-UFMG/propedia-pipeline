rule ml_classifiers:
    input:
        signatures = f"{OUT}/seq_signatures.tsv",
    output:
        report = f"{OUT}/ml_report.tsv",
    params:
        csv_dir      = config["machine"]["ml_train_dir"],
        ifeature_dir = config["machine"]["ifeature_dir"],
        classes      = config["ml"]["classes"],
    script:
        "../scripts/run_ml.py"