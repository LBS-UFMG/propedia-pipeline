rule assemble:
    input:
        pairs     = f"{OUT}/pairs.tsv",
        physchem  = f"{OUT}/physchem.tsv",
        surface   = f"{OUT}/surface.tsv",
        prodigy   = f"{OUT}/prodigy.tsv",
        interface = f"{OUT}/interface.tsv",
        metadata  = f"{OUT}/metadata.tsv",
        clusters    = f"{OUT}/cnr_clusters.tsv",
        therapeutic = f"{OUT}/therapeutic.tsv",
    output:
        csv = f"{OUT}/propedia.csv",
    script:
        "../scripts/assemble.py"