rule cnr_cluster:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        clusters = f"{OUT}/cnr_clusters.tsv",
    script:
        "../scripts/run_cnr.py"

rule legacy_clusters:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        clusters = f"{OUT}/legacy_clusters.tsv",
    params:
        cluster_dir = config["clusters"]["legacy_dir"],
    script:
        "../scripts/legacy_clusters.py"