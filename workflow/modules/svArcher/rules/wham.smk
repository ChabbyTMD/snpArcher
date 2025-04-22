include: "../common.smk"
rule wham_call:
	input:
		unpack(get_bams),
		ref = "results/{refGenome}/data/genome/{refGenome}.fna",
	output:
		wham_vcf=temp("results/{refGenome}/SV/wham/{sample}.vcf"),
	log:
		"logs/{refGenome}/SV/wham/{sample}.log"
	conda:
		"../envs/wham.yaml"
	params:
		contigs = read_contig_file(config["include_contigs"])
	threads: 7
	shell:
		"""
		whamg -x {threads} -c {params.contigs} -a {input.ref} -f {input.bam} > {output.wham_vcf} 2> {log}
		"""
rule wham_vcf_gzip:
	input:
		wham_vcf="results/{refGenome}/SV/wham/{sample}.vcf",
	output:
		wham_vcf_gz = "results/{refGenome}/SV/wham/{sample}.vcf.gz",
		wham_tbi = "results/{refGenome}/SV/wham/{sample}.vcf.gz.tbi",
	conda:
		"../envs/wham.yaml"
	log:
		"logs/{refGenome}/SV/wham/{sample}.vcf.gz.log"
	shell:
		"""
		bgzip -c {input.wham_vcf} > {output.wham_vcf_gz}
		tabix -p vcf {output.wham_vcf_gz}
		"""