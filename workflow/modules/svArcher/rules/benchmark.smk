# Rule to separate DEL DUPS and INV for benchmarking
rule generate_sv_benchmark_test_vcfs:
    input:
        filtered_vcf="results/{refGenome}/SV/postprocess/processed/all_samples_merged.vcf",
    output:
        del_vcf=temp("results/{refGenome}/SV/benchmark/{sample}.del.vcf"),
        dup_vcf=temp("results/{refGenome}/SV/benchmark/{sample}.dup.vcf"),
        inv_vcf=temp("results/{refGenome}/SV/benchmark/{sample}.inv.vcf"),
    conda:
        "../envs/survivor.yaml"
    log:
        "logs/{refGenome}/SV/benchmark/{sample}.log"
    shell:
        """
        bcftools view -i 'SVTYPE==DEL' {input.filtered_vcf} -o {output.del_vcf} 2> {log}
        bcftools view -i 'INFO/SVTYPE="DUP"' {input.filtered_vcf} -o {output.dup_vcf} 2>> {log}
        bcftools view -i 'INFO/SVTYPE="INV"' {input.filtered_vcf} -o {output.inv_vcf} 2>> {log}
        """
rule sort_benchmark_vcfs:
    input:
        del_vcf="results/{refGenome}/SV/benchmark/{sample}.del.vcf",
        dup_vcf="results/{refGenome}/SV/benchmark/{sample}.dup.vcf",
        inv_vcf="results/{refGenome}/SV/benchmark/{sample}.inv.vcf",
    output:
        del_vcf_sorted="results/{refGenome}/SV/benchmark/{sample}.del.sorted.vcf",
        dup_vcf_sorted="results/{refGenome}/SV/benchmark/{sample}.dup.sorted.vcf",
        inv_vcf_sorted="results/{refGenome}/SV/benchmark/{sample}.inv.sorted.vcf",
    conda:
        "../envs/survivor.yaml"
    log:
        "logs/{refGenome}/SV/benchmark/sort_{sample}.log"
    shell:
        """
        bcftools sort {input.del_vcf} -o {output.del_vcf_sorted} 2> {log}
        bcftools sort {input.dup_vcf} -o {output.dup_vcf_sorted} 2>> {log}
        bcftools sort {input.inv_vcf} -o {output.inv_vcf_sorted} 2>> {log}
        """
# Rule to benchamrk SV calls against Arabidopsis Gotkay set
