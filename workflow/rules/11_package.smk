# Package the finished outputs into the propedia26 website file tree.
# Explicit target (not part of `rule all`): `snakemake package`.
# web_dir default = results/<mode>/web (standalone tree to rsync into the site);
# override paths.web_dir to write straight into a propedia26 public/data checkout.
WEB = config["paths"].get("web_dir") or f"{OUT}/web"

rule package:
    input:
        propedia = f"{OUT}/propedia.csv",
        multipro = f"{OUT}/multipro_final.csv",
        mp_cocada = f"{OUT}/multipro_cocada_summary.tsv",
    output:
        marker = f"{WEB}/.packaged",
    params:
        mode                  = config["mode"],
        web_dir               = WEB,
        cocada_dir            = f"{STATE}/cocada",
        multipro_cocada_dir   = f"{STATE}/multipro_cocada",
        legacy_dir            = config["clusters"]["legacy_dir"],
        # optional: a file of target column names -> permutes to the site's order
        column_order          = config.get("package", {}).get("column_order", ""),
        column_order_multipro = config.get("package", {}).get("column_order_multipro", ""),
    script:
        "../scripts/make_package.py"


# Build the website's bulk-download ZIP bundles from the FINISHED outputs.
# Explicit target (not part of `rule all`): `snakemake zips`.
# READ-ONLY over state/ + results/ — creates zips only, never recomputes anything.
rule zips:
    input:
        propedia   = f"{OUT}/propedia.csv",
        multipro   = f"{OUT}/multipro_final.csv",
        seq_sig    = f"{OUT}/seq_signatures.tsv",
        struct_sig = f"{OUT}/struct_signatures.csv",
        pep_marker = f"{STATE}/pep_pdb/.written",
        mp_marker  = f"{STATE}/multipro_cif/.written",
    output:
        marker = f"{WEB}/data/.zipped",
    params:
        out_dir          = f"{WEB}/data",
        cif_dir          = f"{STATE}/cif",
        pep_pdb_dir      = f"{STATE}/pep_pdb",
        multipro_cif_dir = f"{STATE}/multipro_cif",
        legacy_dir       = config["clusters"]["legacy_dir"],
        names            = config.get("package", {}).get("zip_names", {}),
        compresslevel    = config.get("package", {}).get("zip_compresslevel", 6),
    script:
        "../scripts/make_zips.py"
