# Dry run with shell commands printed
snakemake -np --cores 28 --snakefile ../workflow/Snakefile --profile ../profiles/default --conda-frontend conda --use-conda
snakemake --cores 28 --snakefile ../workflow/Snakefile --profile ../profiles/default --conda-frontend conda --use-conda
snakemake -np --cores 28 --snakefile ../workflow/Snakefile --profile ../profiles/default --conda-frontend conda --use-conda --unlock
snakemake --cores 28 --snakefile ../workflow/Snakefile --profile ../profiles/default --conda-frontend conda --use-conda --rerun-incomplete