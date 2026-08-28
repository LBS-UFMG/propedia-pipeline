rule physchem:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        physchem = f"{OUT}/physchem.tsv",
    script:
        "../scripts/physchem.py"

rule prodigy:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        prodigy = f"{OUT}/prodigy.tsv",
        errors  = f"{OUT}/prodigy_errors.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir         = f"{STATE}/cif",
        distance_cutoff = config["prodigy"]["distance_cutoff"],
        temperature     = config["prodigy"]["temperature"],
        ckpt            = f"{STATE}/checkpoint/prodigy",
    script:
        "../scripts/run_prodigy.py"

rule cocada:
    input:
        pairs = f"{OUT}/pairs.tsv",
        done  = f"{OUT}/interim/download_cifs.done",
    output:
        summary = f"{OUT}/cocada_summary.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir    = config["machine"]["cif_dir"],
        cocada_dir = config["machine"]["cocada_dir"],
        ph         = config["cocada"]["ph"],
        out_root   = f"{STATE}/cocada",
        ckpt       = f"{STATE}/checkpoint/cocada",
    script:
        "../scripts/run_cocada.py"

rule ifeature:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        signatures = f"{OUT}/seq_signatures.tsv",
        excluded   = f"{OUT}/seq_signatures_excluded.tsv",
    threads: config["compute"]["threads"]
    params:
        ifeature_dir      = config["machine"]["ifeature_dir"],
        include_ctriad    = False,
        min_signature_len = config["signatures"]["min_signature_len"],
        ckpt              = f"{STATE}/checkpoint/ifeature",
    script:
        "../scripts/run_ifeature.py"

rule peptide_pdbs:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        marker = f"{STATE}/pep_pdb/.written",
    params:
        pair_cif_dir = f"{STATE}/cif",
        pep_pdb_dir  = f"{STATE}/pep_pdb",
    script:
        "../scripts/write_peptide_pdbs.py"

rule signa:
    input:
        marker = f"{STATE}/pep_pdb/.written",
    output:
        signatures = f"{OUT}/struct_signatures.csv",
    threads: config["compute"]["threads"]
    params:
        signa_dir    = config["machine"]["signa_dir"],
        pep_pdb_dir  = f"{STATE}/pep_pdb",
        ckpt         = f"{STATE}/checkpoint/signa",
        cutoff_limit = config["signatures"]["acsm_cutoff_limit"],
        cutoff_step  = config["signatures"]["acsm_cutoff_step"],
        cumulative   = config["signatures"]["acsm_cumulative"],
    script:
        "../scripts/run_signa.py"

rule freesasa:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        surface = f"{OUT}/surface.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir       = f"{STATE}/cif",
        bsa_threshold = config["cutoffs"]["bsa_threshold_A2"],
        ckpt          = f"{STATE}/checkpoint/freesasa",
    script:
        "../scripts/run_freesasa.py"

rule interface:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        interface = f"{OUT}/interface.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir = f"{STATE}/cif",
        cutoff  = config["cutoffs"]["interaction_distance_A"],
        ckpt    = f"{STATE}/checkpoint/interface",
    script:
        "../scripts/interface_residues.py"

rule metadata:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        metadata = f"{OUT}/metadata.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir = config["machine"]["cif_dir"],
        ckpt    = f"{STATE}/checkpoint/metadata",
    script:
        "../scripts/metadata.py"