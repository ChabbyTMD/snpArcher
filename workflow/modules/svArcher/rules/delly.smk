
rule delly_call:
    input:
        unpack(get_bams),
        ref = "results/{refGenome}/data/genome/{refGenome}.fna",
    output:
        delly_bcf=temp("results/{refGenome}/SV/delly/{sample}.bcf"),
        delly_csi=temp("results/{refGenome}/SV/delly/{sample}.bcf.csi"),
    log:
        "logs/{refGenome}/SV/delly/{sample}.log"
    conda:
        "../envs/delly.yaml"
    shell:
        "delly call -g {input.ref} -o {output.delly_bcf} {input.bam} 2> {log}"


rule delly_vcf:
    input:
        delly_bcf="results/{refGenome}/SV/delly/{sample}.bcf",
    output:
        delly_vcf="results/{refGenome}/SV/delly/{sample}.vcf.gz",
        delly_tbi="results/{refGenome}/SV/delly/{sample}.vcf.gz.tbi",
    conda:
        "../envs/delly.yaml"
    log:
        "logs/{refGenome}/SV/delly/{sample}_vcf.log"
    shell:
        """
        bcftools view {input} -O z -o {output.delly_vcf} 2> {log}
        tabix -p vcf {output.delly_vcf}
        """