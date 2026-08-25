rule validate:
    input:
        pairs      = f"{OUT}/pairs.tsv",
        physchem   = f"{OUT}/physchem.tsv",
        prodigy    = f"{OUT}/prodigy.tsv",
        cocada     = f"{OUT}/cocada_summary.tsv",
        seq_sig    = f"{OUT}/seq_signatures.tsv",
        surface    = f"{OUT}/surface.tsv",        
        struct_sig = f"{OUT}/struct_signatures.csv",
        clusters   = f"{OUT}/cnr_clusters.tsv",
        ml         = f"{OUT}/ml_report.tsv",
        
    output:
        report = f"{OUT}/reproduction_report.txt",
    params:
        oracle_csv = config["machine"]["oracle_csv"],
        mode       = config["mode"],
    script:
        "../scripts/validate.py"
        