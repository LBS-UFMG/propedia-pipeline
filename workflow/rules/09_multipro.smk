rule multipro_extract:
    input:
        pairs   = f"{OUT}/pairs.tsv",
        surface = f"{OUT}/surface.tsv",
    output:
        multipro = f"{OUT}/multipro.tsv",
    script:
        "../scripts/multipro_extract.py"
		
rule multipro_cifs:
    input:
        multipro = f"{OUT}/multipro.tsv",
    output:
        marker = f"{OUT}/multipro_cif/.written",
    params:
        cif_dir     = config["machine"]["cif_dir"],
        cif_out_dir = f"{OUT}/multipro_cif",
    script:
        "../scripts/multipro_cifs.py"

rule multipro_surface:
    input:
        multipro = f"{OUT}/multipro.tsv",
        marker   = f"{OUT}/multipro_cif/.written",
    output:
        surface = f"{OUT}/multipro_surface.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir = f"{OUT}/multipro_cif",
        ckpt    = f"{OUT}/.checkpoint/multipro_surface",
    script:
        "../scripts/multipro_surface.py"

rule multipro_prodigy:
    input:
        multipro = f"{OUT}/multipro.tsv",
        marker   = f"{OUT}/multipro_cif/.written",
    output:
        prodigy = f"{OUT}/multipro_prodigy.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir         = f"{OUT}/multipro_cif",
        distance_cutoff = config["prodigy"]["distance_cutoff"],
        temperature     = config["prodigy"]["temperature"],
        ckpt            = f"{OUT}/.checkpoint/multipro_prodigy",
    script:
        "../scripts/multipro_prodigy.py"

rule multipro_assemble:
    input:
        multipro  = f"{OUT}/multipro.tsv",
        surface   = f"{OUT}/multipro_surface.tsv",
        prodigy   = f"{OUT}/multipro_prodigy.tsv",
        physchem  = f"{OUT}/physchem.tsv",
        metadata  = f"{OUT}/metadata.tsv",
        interface = f"{OUT}/interface.tsv",
        legacy    = f"{OUT}/legacy_clusters.tsv",
    output:
        csv = f"{OUT}/multipro_final.csv",
    script:
        "../scripts/multipro_assemble.py"		