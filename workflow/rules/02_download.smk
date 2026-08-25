rule download_cifs:
    input:
        ids = f"{OUT}/interim/candidate_ids.txt",
    output:
        done     = f"{OUT}/interim/download_cifs.done",
        failures = f"{OUT}/interim/download_failures.txt",
    threads: config["download"]["workers"]
    params:
        cif_dir         = config["machine"]["cif_dir"],
        max_retries     = config["download"]["max_retries"],
        timeout_seconds = config["download"]["timeout_seconds"],
    script:
        "../scripts/download_cifs.py"