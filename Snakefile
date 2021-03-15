AGS = ["033",
      "034",
      "035",
      "036",
      "037",
      "038",
      "039",
      "040",
      "041",
      "042",
      "043",
      "044",
      "045",
      "046",
      "047",
      "048",
      "049",
      "050"]
ACS = ["001",
      "002",
      "003",
      "006",
      "005",
      "007",
      "008",
      "009",
      "010",
      "011",
      "012",
      "081",
      "082",
      "083",
      "084",
      "085",
      "086",
      "087",
      "088",
      "089",
      "090"]

samplesALL = AGS + ACS
bwa_path = "bwa"
samtools_path = "samtools"
gatk_path = "gatk"
picard_path = "picard"
bedtools_path = "bedtools"

rule all:
	input:
		"howlerbams/howlerbams/gatk.called.raw.vcf.gz"

rule revert_bam_AG:
	input:
		bamAG = "howlerbams/IonXpress_{AG}_R_2018_08_09_13_40_32_user_Eve-167-2018_08_09_Nicole_2_Auto_user_Eve-167-2018_08_09_Nicole_2_225.bam"
	output:
		unaligned_bam = "howlerbams/howlerbams/{AG}.un.bam"
	params:
		picard = picard_path
	shell:
		"{params.picard} RevertSam I={input.bamAG} O={output.unaligned_bam}"

rule revert_bam_AC:
	input:
		bamAC = "howlerbams/IonXpress_{AC}_R_2018_04_18_12_41_37_user_Eve-162-2018_04_18_Nicole_1_Auto_user_Eve-162-2018_04_18_Nicole_1_220.bam"
	output:
		unaligned_bam = "howlerbams/howlerbams/{AC}.un.bam"
	params:
		picard = picard_path
	shell:
		"{params.picard} RevertSam I={input.bamAC} O={output.unaligned_bam}"

rule bam_to_fastq:
	input:
		unaligned_bam = "howlerbams/{samples}.un.bam"
	output:
		fastq = "howlerbams/howlerbams/{samples}.fq"
	params:
		bedtools = bedtools_path
	shell:
		"{params.bedtools} bamtofastq -i {input.unaligned_bam} -fq {output.fastq}"


rule map_fq:
	input:
		fq =  "howlerbams/howlerbams/{samples}.fq",
		ref  = "howlerbams/ACrefs/ACref97.fasta"
	output:
		bams = "howlerbams/howlerbams/{samples}.sorted.bam"
	params:
		id = "{samples}",
		pl = "IonPGM",
		bwa = bwa_path,
		samtools = samtools_path
	shell:
		"{params.bwa} mem -R "
		"'@RG\\tID:{params.id}\\tSM:{params.id}\\tLB:{params.id}\\tPU:{params.id}\\tPL:{params.pl}' "
		"{input.ref} {input.fq} "
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output.bams}"

rule index_sorted_bam:
	input:
		"howlerbams/howlerbams/{samples}.sorted.bam"
	output:
		"howlerbams/howlerbams/{samples}.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "howlerbams/howlerbams/{samples}.sorted.bam",
		bai = "howlerbams/howlerbams/{samples}.sorted.bam.bai"
	output:
		bam = "howlerbams/howlerbams/{samples}.sorted.mkdup.bam",
		metrics = "howlerbams/howlerbams/{samples}.picard_mkdup_metrics.txt"
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx1g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam:
	input:
		"howlerbams/howlerbams/{samples}.sorted.mkdup.bam"
	output:
		"howlerbams/howlerbams/{samples}.sorted.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule gatk_gvcf:
	input:
		ref = "howlerbams/ACrefs/ACref97.fasta",
		bam = "howlerbams/howlerbams/{samples}.sorted.mkdup.bam",
		bai = "howlerbams/howlerbams/{samples}.sorted.mkdup.bam.bai"
	output:
		"howlerbams/howlerbams/{samples}.g.vcf.gz"
	params:
		gatk = gatk_path
	shell:
		"""{params.gatk} --java-options "-Xmx1g" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} """
		"""-ERC GVCF -O {output}"""

rule gatk_combinegvcfs:
	input:
		ref = "howlerbams/ACrefs/ACref97.fasta",
		gvcfs = lambda wildcards: expand(
			"howlerbams/howlerbams/{samples}.g.vcf.gz",
			sample=samplesALL)
	output:
		"howlerbams/howlerbams/gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		print(
			"""{params.gatk} --java-options "-Xmx1g" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}"""
		)
		shell(
			"""{params.gatk} --java-options "-Xmx1g" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}""")

rule gatk_genotypegvcf:
	input:
		ref = "howlerbams/ACrefs/ACref97.fasta",
		gvcf = "howlerbams/howlerbams/gatk.combinegvcf.g.vcf.gz"
	output:
		"howlerbams/howlerbams/gatk.called.raw.vcf.gz"
	params:
		gatk = gatk_path
	shell:
		"""{params.gatk} --java-options "-Xmx1g" """
		"""GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"""
