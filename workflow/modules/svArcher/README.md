# *svArcher*: A Structural Variant Calling Extension module for snpArcher.

> Under construction*.


*svArcher* is a reproducible subworkflow designed to call structural variants (SVs) from short read genomic sequencing data. The pipeline utilizes an ensemble of structural variant callers (DELLY, Lumpy and Wham) that call SV's independently. Results are merged and filtered with SURVIVOR. Current parameters (500 1 1 0 50) are set to achieve a concensus of SV's from all 3 callers.

## Usage

svArcher is deactivated by default. Set svCall option to `True` in your workflow `config.yaml` file.

Additionally, provide a `.csv` file with a list of contigs/chromosome ID's present in your reference file

```yaml
svCall: False # Set to true to activate svArcher
include_contigs: "/path/to/ref_contigs.csv"
```

