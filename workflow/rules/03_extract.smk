rule build_sample:
    input:
        oracle = config["machine"]["oracle_csv"],
        done   = f"{OUT}/interim/download_cifs.done",
    output:
        sample = f"{OUT}/interim/sample_ids.txt",
    params:
        cif_dir     = config["machine"]["cif_dir"],
        sample_size = config["sample"]["size"],
        seed        = config["sample"]["seed"],
    script:
        "../scripts/build_sample.py"

def extract_id_source(wildcards):
    if config["mode"] == "sample":
        return f"{OUT}/interim/sample_ids.txt"
    return f"{OUT}/interim/candidate_ids.txt"

rule extract_pairs:
    input:
        ids  = extract_id_source,
        done = f"{OUT}/interim/download_cifs.done",
    output:
        pairs = f"{OUT}/pairs.tsv",
    threads: config["compute"]["threads"]
    params:
        cif_dir     = config["machine"]["cif_dir"],
        cif_out_dir = f"{STATE}/cif",
        pep_min     = config["cutoffs"]["peptide_len_min"],
        pep_max     = config["cutoffs"]["peptide_len_max"],
        prot_min    = config["cutoffs"]["protein_len_min"],
        cutoff      = config["cutoffs"]["interaction_distance_A"],
        max_atoms   = config["cutoffs"].get("max_atoms_per_structure", 0),
        atom_cache  = f"{STATE}/atom_counts.tsv",
        oversized   = f"{OUT}/oversized.tsv",
        ckpt        = f"{STATE}/checkpoint/extract",
    script:
        "../scripts/extract_pairs.py"