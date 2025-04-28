rule discordant_extract:
    input:
        unpack(get_bams),
    output:
        unsorted_discordant_bam = temp("results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.unsorted.discordant.bam"),
    log:
        "logs/{refGenome}/SV/lumpy_prep/{sample}.discordant_extract.txt"
    conda:
        "../envs/lumpy.yaml"
    shell:
        """
        samtools view -b -F 1294 {input.bam} > {output.unsorted_discordant_bam} 2> {log}
        """

rule split_read_extract:
    input:
        unpack(get_bams),
    output:
        unsorted_split_read_bam = temp("results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.unsorted.split_read.bam"),
    log:
        "logs/{refGenome}/SV/lumpy_prep/{sample}.split_read_extract.txt"
    conda:
        "../envs/lumpy.yaml"
    shell:
        """
        samtools view -h {input.bam} \
            | extractSplitReads_BwaMem -i stdin \
            | samtools view -Sb - \
            > {output.unsorted_split_read_bam} 2> {log}
        """

rule sort_splits_discordants:
    input:
        unsorted_discordant_bam = "results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.unsorted.discordant.bam",
        unsorted_split_read_bam = "results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.unsorted.split_read.bam",
    output:
        discordant_bam = temp("results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.discordant.bam"),
        split_read_bam = temp("results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.split_read.bam"),
    log:
        "logs/{refGenome}/SV/lumpy_prep/{sample}.sort_splits_discordants.txt"
    conda:
        "../envs/lumpy.yaml"
    shell:
        """
        samtools sort -o {output.discordant_bam} {input.unsorted_discordant_bam} 2> {log}
        samtools sort -o {output.split_read_bam} {input.unsorted_split_read_bam} 2>> {log}
        """

rule lumpy_call:
    """
    Rule that implements the split read methodology using lumpy
    """
    input:
        unpack(get_bams),
        fai = "results/{refGenome}/data/genome/{refGenome}.fna.fai",
        sorted_discordant_bam = "results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.discordant.bam",
        sorted_split_read_bam = "results/{refGenome}/SV/lumpy_prep/{sample}/{sample}.split_read.bam",
    output:
        lumpy_vcf = temp("results/{refGenome}/SV/lumpy/{sample}.vcf"),
        lumpy_vcf_gz = "results/{refGenome}/SV/lumpy/{sample}.vcf.gz",
    log:
        "logs/{refGenome}/SV/lumpy/{sample}.txt"
    conda:
        "../envs/lumpy.yaml"
    shell:
    # Call SVs with lumpy express, compress, reheader and index SV call file # changed
    # NOTE: check vcf file after calling whether chromosome names are correct
        """
        lumpyexpress \
            -B {input.bam} \
            -S {input.sorted_split_read_bam} \
            -D {input.sorted_discordant_bam} \
            -o {output.lumpy_vcf} 2> {log}
        bgzip -c {output.lumpy_vcf} > {output.lumpy_vcf_gz}
        """