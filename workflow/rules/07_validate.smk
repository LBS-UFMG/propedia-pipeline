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
        propedia   = f"{OUT}/propedia.csv",
        multipro   = f"{OUT}/multipro_final.csv",
    output:
        report = f"{OUT}/reproduction_report.txt",
    params:
        oracle_csv      = config["machine"]["oracle_csv"],
        mode            = config["mode"],
        multipro_oracle = config["machine"]["multipro_oracle_csv"],
    script:
        "../scripts/validate.py"
