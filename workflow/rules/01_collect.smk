rule fetch_candidate_ids:
    output:
        ids     = f"{OUT}/interim/candidate_ids.txt",
        meta    = f"{OUT}/interim/candidate_ids.meta.json",
        missing = f"{OUT}/interim/candidate_ids.oracle_missing.txt",
    params:
        snapshot_date        = config["pdb_snapshot_date"],
        protein_polymer_type = config["collection"]["protein_polymer_type"],
        protein_len_min      = config["cutoffs"]["protein_len_min"],
        min_polymer_chains   = config["collection"]["min_polymer_chains"],
        oracle_ids           = config["machine"]["oracle_csv"],
    script:
        "../scripts/fetch_pdb_ids.py"