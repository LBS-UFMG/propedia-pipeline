rule provenance:
    input:
        surface = f"{OUT}/surface.tsv",
    output:
        provenance = f"{OUT}/provenance.tsv",
    params:
        snapshot_date = config["pdb_snapshot_date"],
        ledger        = f"{STATE}/provenance.tsv",     # DURABLE: preserve between releases
    script:
        "../scripts/provenance.py"

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
        legacy      = f"{OUT}/legacy_clusters.tsv",
        provenance  = f"{OUT}/provenance.tsv",
        pisa        = f"{OUT}/pisa.tsv",
    output:
        csv = f"{OUT}/propedia.csv",
    script:
        "../scripts/assemble.py"