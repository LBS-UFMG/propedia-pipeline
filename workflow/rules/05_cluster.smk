rule cnr_cluster:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        clusters = f"{OUT}/cnr_clusters.tsv",
    script:
        "../scripts/run_cnr.py"