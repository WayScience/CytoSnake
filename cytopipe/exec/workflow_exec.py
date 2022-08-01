from pathlib import Path
import snakemake

def exec_preprocessing(params_list: list) -> int:

    # get workflow path
    snakefile = str(Path("Snakefile").absolute())

    # execute
    status = snakemake.snakemake(snakefile, use_conda=True)

    if status:
        return 0
    return 1
