# Snapshot the finished deliverables into a dated, immutable release directory with a
# provenance manifest. Explicit target (not part of `rule all`): `snakemake release`.
rule release:
    input:
        propedia = f"{OUT}/propedia.csv",
        multipro = f"{OUT}/multipro_final.csv",
        report   = f"{OUT}/reproduction_report.txt",
    output:
        manifest = f"releases/propedia-{config['pdb_snapshot_date']}/manifest.json",
    params:
        release_dir   = f"releases/propedia-{config['pdb_snapshot_date']}",
        snapshot_date = config["pdb_snapshot_date"],
        mode          = config["mode"],
        state_dir     = STATE,
    script:
        "../scripts/make_release.py"
