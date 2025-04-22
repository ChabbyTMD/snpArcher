import sys
from pathlib import Path
import pandas as pd
import snparcher_utils
from pkg_resources import parse_version
# Get utils. This is not great, but we can move to setup.py and install later if want
utils_path = (Path(workflow.main_snakefile).parent.parent.parent).resolve()
if str(utils_path) not in sys.path:
    sys.path.append(str(utils_path))

# Can't be less than 7 cuz of min version in snakefile
SNAKEMAKE_VERSION = 8 if parse_version(snakemake.__version__) >= parse_version("8.0.0") else 7
logger.warning(f"snpArcher: Using Snakemake {snakemake.__version__}")
if SNAKEMAKE_VERSION >= 8:
    DEFAULT_STORAGE_PREFIX = StorageSettings.default_storage_prefix if StorageSettings.default_storage_prefix is not None else ""
else:
    # backwards compatability w/ snakemake <= 7
    DEFAULT_STORAGE_PREFIX = workflow.default_remote_prefix
    if config["remote_reads"]:
        from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
        GS = GSRemoteProvider()
        GS_READS_PREFIX = config['remote_reads_prefix']


def get_bams(wc):
    out = {"bam": None, "bai": None}
    if config["mark_duplicates"]:
        out["bam"] = "results/{refGenome}/bams/{sample}_final.bam"
        out["bai"] = "results/{refGenome}/bams/{sample}_final.bam.bai"
        return out
    else:
        return dedup_input(wc)

def read_contig_file(contig_file: str) -> list[str]:
    """
    Reads a one column csv file containing contig names a user would want to call SV's into a list of contig names.
    
    """
    with open(contig_file, "r") as file:
        contigs = [line.strip() for line in file]
    return contigs

# TODO: Define function to handle svArcher outputs.

def svArcher_output(wildcards):
    output = []
    output.extend(expand("results/{refGenome}/SV/{method}/{sample}.vcf.gz", refGenome=REFGENOME, method=["delly", "lumpy", "wham"], sample=samples["BioSample"].unique().tolist()))
    return output