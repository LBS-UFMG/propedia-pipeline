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
    params:
        pdb_dir         = f"{OUT}/pdb",
        distance_cutoff = config["prodigy"]["distance_cutoff"],
        temperature     = config["prodigy"]["temperature"],
    script:
        "../scripts/run_prodigy.py"

rule cocada:
    input:
        pairs = f"{OUT}/pairs.tsv",
        done  = f"{OUT}/interim/download_cifs.done",
    output:
        summary = f"{OUT}/cocada_summary.tsv",
    params:
        cif_dir    = config["machine"]["cif_dir"],
        cocada_dir = config["machine"]["cocada_dir"],
        ph         = config["cocada"]["ph"],
        out_root   = f"{OUT}/cocada",
    script:
        "../scripts/run_cocada.py"

rule ifeature:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        signatures = f"{OUT}/seq_signatures.tsv",
        excluded   = f"{OUT}/seq_signatures_excluded.tsv",
    params:
        ifeature_dir      = config["machine"]["ifeature_dir"],
        include_ctriad    = False,
        min_signature_len = config["signatures"]["min_signature_len"],
    script:
        "../scripts/run_ifeature.py"

rule peptide_pdbs:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        marker = f"{OUT}/pep_pdb/.written",
    params:
        pair_pdb_dir = f"{OUT}/pdb",
        pep_pdb_dir  = f"{OUT}/pep_pdb",
    script:
        "../scripts/write_peptide_pdbs.py"

rule signa:
    input:
        marker = f"{OUT}/pep_pdb/.written",
    output:
        signatures = f"{OUT}/struct_signatures.csv",
    params:
        signa_dir    = config["machine"]["signa_dir"],
        pep_pdb_dir  = f"{OUT}/pep_pdb",
        raw_output   = f"{OUT}/signa_raw.csv",
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
    params:
        pdb_dir       = f"{OUT}/pdb",
        bsa_threshold = config["cutoffs"]["bsa_threshold_A2"],
    script:
        "../scripts/run_freesasa.py"
		
rule interface:
    input:
        pairs = f"{OUT}/pairs.tsv",
    output:
        interface = f"{OUT}/interface.tsv",
    params:
        pdb_dir = f"{OUT}/pdb",
        cutoff  = config["cutoffs"]["interaction_distance_A"],
    script:
        "../scripts/interface_residues.py"		