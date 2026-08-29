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
        marker = f"{STATE}/multipro_cif/.written",
    params:
        cif_dir     = config["machine"]["cif_dir"],
        cif_out_dir = f"{STATE}/multipro_cif",
    script:
        "../scripts/multipro_cifs.py"

rule multipro_surface:
    input:
        multipro = f"{OUT}/multipro.tsv",
        marker   = f"{STATE}/multipro_cif/.written",
    output:
        surface = f"{OUT}/multipro_surface.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir = f"{STATE}/multipro_cif",
        ckpt    = f"{STATE}/checkpoint/multipro_surface",
    script:
        "../scripts/multipro_surface.py"

rule multipro_cocada:
    input:
        multipro = f"{OUT}/multipro.tsv",
        done     = f"{OUT}/interim/download_cifs.done",   # runs on the ORIGINAL CIFs
    output:
        summary = f"{OUT}/multipro_cocada_summary.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir    = config["machine"]["cif_dir"],        # original RCSB CIFs (COCaDA can't
                                                          # parse the MMCIFIO multipro CIF)
        cocada_dir = config["machine"]["cocada_dir"],
        ph         = config["cocada"]["ph"],
        out_root   = f"{STATE}/multipro_cocada",
        ckpt       = f"{STATE}/checkpoint/multipro_cocada",
    script:
        "../scripts/multipro_cocada.py"

rule multipro_prodigy:
    input:
        multipro = f"{OUT}/multipro.tsv",
        marker   = f"{STATE}/multipro_cif/.written",
    output:
        prodigy = f"{OUT}/multipro_prodigy.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir         = f"{STATE}/multipro_cif",
        distance_cutoff = config["prodigy"]["distance_cutoff"],
        temperature     = config["prodigy"]["temperature"],
        ckpt            = f"{STATE}/checkpoint/multipro_prodigy",
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
        provenance = f"{OUT}/provenance.tsv",
        pisa      = f"{OUT}/pisa.tsv",            # pep-pro PISA, aggregated per multipro entry
    output:
        csv = f"{OUT}/multipro_final.csv",
    script:
        "../scripts/multipro_assemble.py"		