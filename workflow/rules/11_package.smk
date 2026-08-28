# Package the finished outputs into the propedia26 website file tree.
# Explicit target (not part of `rule all`): `snakemake package`.
# web_dir default = results/<mode>/web (standalone tree to rsync into the site);
# override paths.web_dir to write straight into a propedia26 public/data checkout.
WEB = config["paths"].get("web_dir") or f"{OUT}/web"

rule package:
    input:
        propedia = f"{OUT}/propedia.csv",
        multipro = f"{OUT}/multipro_final.csv",
    output:
        marker = f"{WEB}/.packaged",
    params:
        mode                  = config["mode"],
        web_dir               = WEB,
        cocada_dir            = f"{STATE}/cocada",
        legacy_dir            = config["clusters"]["legacy_dir"],
        # optional: a file of target column names -> permutes to the site's order
        column_order          = config.get("package", {}).get("column_order", ""),
        column_order_multipro = config.get("package", {}).get("column_order_multipro", ""),
    script:
        "../scripts/make_package.py"
