# *svArcher*: A Structural Variant Calling Extension module for snpArcher.

<div align="center">
<img src="img/svArcher_Logo.png" alt="svArcher Logo" width="150" height="190">
</div>

*svArcher* is a reproducible subworkflow designed to call structural variants (SVs) from short read genomic sequencing data. The pipeline utilizes an ensemble of structural variant callers (DELLY, Lumpy and Wham) that call SV's independently. Results are merged and filtered with SURVIVOR. Current parameters (1000 3 1 1 0 50) are set to achieve a concensus of SV's from all 3 callers per sample.

## Usage

svArcher is deactivated by default. Set svCall option to `True` in your workflow `config.yaml` file.

Additionally, provide a `.csv` file with a list of contigs/chromosome ID's present in the reference sequence file. Each ID must start on its own line.

```yaml
svCall: False # Set to true to activate svArcher
include_contigs: "/path/to/ref_contigs.csv"
```

## Benchmarking

This module enables users to generate benchmark comparison metrics for deletion, duplication and inversion structural variant types. To activate this feature, set svBenchmark to `True`.

We recommend creating a single directory to store all ground truth `.vcf` files. When configuring options below, provide the absolute path to this directory.

```yaml
svBenchmark: False
# Benchmark Truth set paths
deletions: "benchmark_truth/"
duplications: "benchmark_truth/"
inversions: "benchmark_truth/"

```
