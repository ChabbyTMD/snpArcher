# Sort vcf files from delly, lumpy and wham

rule sort_vcfs:
	input:
		delly_vcf="results/{refGenome}/SV/delly/{sample}.vcf",
		lumpy_vcf="results/{refGenome}/SV/lumpy/{sample}.vcf",
		wham_vcf="results/{refGenome}/SV/wham/{sample}.vcf",
	output:
		delly_vcf_sorted="results/{refGenome}/SV/delly/{sample}.sorted.vcf",
		lumpy_vcf_sorted="results/{refGenome}/SV/lumpy/{sample}.sorted.vcf",
		wham_vcf_sorted="results/{refGenome}/SV/wham/{sample}.sorted.vcf",
	conda:
		"../envs/wham.yaml"
	log:
		"logs/{refGenome}/SV/postprocess/sort_{sample}.log"
	shell:
		"""
		bcftools sort {input.delly_vcf} -o {output.delly_vcf_sorted} 2> {log}
		bcftools sort {input.lumpy_vcf} -o {output.lumpy_vcf_sorted} 2>> {log}
		bcftools sort {input.wham_vcf} -o {output.wham_vcf_sorted} 2>> {log}
		"""
# Rule to merge SVS from all methods
rule sample_sv_call_merge:
	input:
		delly_vcf="results/{refGenome}/SV/delly/{sample}.sorted.vcf",
		lumpy_vcf="results/{refGenome}/SV/lumpy/{sample}.sorted.vcf",
		wham_vcf="results/{refGenome}/SV/wham/{sample}.sorted.vcf",
	output:
		merged_vcf="results/{refGenome}/SV/postprocess/raw_merge/{sample}.vcf",
	conda:
		"../envs/survivor.yaml"
	log:
		"logs/{refGenome}/SV/postprocess/{sample}.log"
	params:
		sample_file = "results/{refGenome}/SV/postprocess/raw_merge/{sample}.txt",
	shell:
		"""
		# Create a file with the paths to the vcf files
		# This is needed for survivor merge
		ls {input.delly_vcf} {input.lumpy_vcf} {input.wham_vcf} > {params.sample_file}
		SURVIVOR merge {params.sample_file} 1000 3 1 1 0 50 {output.merged_vcf} 2> {log}
		"""
# Rule to filter SV calls from all methods
rule filter_sv_calls:
    input:
        merged_vcf="results/{refGenome}/SV/postprocess/raw_merge/{sample}.vcf",
    output:
        filtered_vcf=temp("results/{refGenome}/SV/postprocess/filtered/{sample}.filtered.vcf"),

    conda:
        "../envs/survivor.yaml"
    log:
        "logs/{refGenome}/SV/postprocess/{sample}.filter.log"
    shell:
        """
		# Filter out SV calls with length > 10kb
        bcftools view -i 'SVLEN <= 10000' {input.merged_vcf} -o {output.filtered_vcf} 2> {log}
        """
rule post_filter_process:
	input:
		filtered_vcf="results/{refGenome}/SV/postprocess/filtered/{sample}.filtered.vcf",
	output:
		processed_vcf="results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf",
	conda:
		"../envs/survivor.yaml"
	log:
		"logs/{refGenome}/SV/postprocess/{sample}.post_filter.log"
	shell:
		"""
		awk '
		BEGIN {{ FS = OFS = "\t" }}
		# Print header lines unchanged
		/^#/ {{ print; next }}
		# Process records where the 5th field is NOT <DEL>, <DUP>, or <INV>
		$5 != "<DEL>" && $5 != "<DUP>" && $5 != "<INV>" {{
			# Find SVTYPE=value in the INFO field ($8)
			if ($8 ~ /SVTYPE=[^;]+/) {{
				# Extract the SVTYPE value
				split($8, info_fields, ";")
				for (i in info_fields) {{
					if (info_fields[i] ~ /^SVTYPE=/) {{
						svtype_value = substr(info_fields[i], 8)
						break
					}}
				}}
				# Replace ALT field ($5) with the extracted SVTYPE value enclosed in <>
				$5 = "<" svtype_value ">"
			}}
			# Print the modified line
			print
			next
		}}
		# Print other records unchanged
		{{ print }}
		' {input.filtered_vcf} > {output.processed_vcf} 2> {log}
		"""

# Define a default value for refGenome if not provided
REF_GENOME = config.get("refGenome", "GCF_000001735.4")

# Update rules to use the default value
rule all_sample_merge:
    input:
        processed_vcfs=expand("results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf", refGenome=REF_GENOME, sample=glob_wildcards("results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf").sample),
    output:
        "results/{refGenome}/SV/postprocess/processed/all_samples_merged.vcf",
    conda:
        "../envs/survivor.yaml",
    log:
        "logs/{refGenome}/SV/postprocess/all_samples.log",
    params:
        sample_vcfs = "results/{refGenome}/SV/postprocess/processed/all_samples.txt",
    shell:
        """
        # Merge all processed vcf files into one
        ls {input.processed_vcfs} > {params.sample_vcfs}
        SURVIVOR merge {params.sample_vcfs} 1000 1 1 1 0 50 {output} 2> {log}
        """
# 
 # processed_vcfs=expand("results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf", refGenome=REF_GENOME, sample=glob_wildcards("results/{refGenome}/SV/postprocess/processed/{sample}.processed.vcf").sample),