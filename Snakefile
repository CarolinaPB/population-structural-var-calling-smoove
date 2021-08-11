configfile: "config.yaml"

from snakemake.utils import makedirs


#################################
# author: Carolina Pita Barros  #
# carolina.pitabarros@wur.nl    #
# date: August 2021             #
#################################

if "OUTDIR" in config:
    print("\nSaving to " + config["OUTDIR"] + "\n")
    workdir: config["OUTDIR"]

makedirs("logs_slurm")

pipeline = "population-structural-var-calling-smoove" # replace with your pipeline's name


include: "rules/create_file_log.smk"

REFERENCE=config["REFERENCE"]
READS_DIR=config["READS_DIR"]
PREFIX = config["PREFIX"]
SAMPLES_LIST = config["SAMPLE_LIST"]
CONTIGS_IGNORE = config["CONTIGS_IGNORE"]

with open(SAMPLES_LIST, "r") as infile:
    content = infile.readlines()
    SAMPLES = [os.path.splitext(x)[0] for x in content]

localrules: create_file_log

rule all:
    input:
        files_log,
        expand("3_genotyped/{sample}-smoove.genotyped.vcf.gz", sample=SAMPLES),
        expand("4_paste/{prefix}.smoove.square.vcf.gz", prefix = PREFIX)



with open(CONTIGS_IGNORE, "r") as infile:
    content =  infile.read().splitlines()
    CONTIGS =  ",".join(content)

rule smoove_call:
    input:
        bam = os.path.join(READS_DIR, "{sample}.bam"),
        reference = REFERENCE,
    output:
        temp("1_call/{sample}-smoove.genotyped.vcf.gz")
    message:
        'Rule {rule} processing'
    params:
        outdir = "1_call",
        contigs=CONTIGS
    conda:
        "envs/smoove.yaml"
    shell:
        """
smoove call --outdir {params.outdir} \
--name {wildcards.sample} \
--fasta {input.reference} \
--duphold \
--genotype \
-p 1 \
--excludechroms {params.contigs} \
{input.bam}
        """

rule smoove_merge:
    input:
        reference = REFERENCE,
        vcfs = expand("1_call/{sample}-smoove.genotyped.vcf.gz", sample=SAMPLES),
    output:
        temp(expand("2_merged/{prefix}.sites.vcf.gz", prefix=PREFIX))
    message:
        'Rule {rule} processing'
    params:
        outdir = "2_merged",
        name = PREFIX
    conda:
        "envs/smoove.yaml"
    shell:
        """
smoove merge --name {params.name} \
--outdir {params.outdir} \
-f {input.reference} \
{input.vcfs}
        """

rule smoove_genotype:
    input:
        reference = REFERENCE,
        vcf = rules.smoove_merge.output,
        bam = os.path.join(READS_DIR, "{sample}.bam"),
    output:
        temp("3_genotyped/{sample}-smoove.genotyped.vcf.gz")
    message:
        'Rule {rule} processing'
    params:
        outdir = "3_genotyped",
        name = "{sample}"
    conda:
        "envs/smoove.yaml"
    shell:
        """
smoove genotype -x -p 3 \
--name {params.name} \
--outdir {params.outdir} \
--fasta {input.reference} \
--duphold \
--vcf {input.vcf} \
{input.bam}
        """

rule smoove_paste:
    input:
        expand("3_genotyped/{sample}-smoove.genotyped.vcf.gz", prefix=PREFIX, sample=SAMPLES),
    output:
        "4_paste/{prefix}.smoove.square.vcf.gz"
    message:
        'Rule {rule} processing'
    conda:
        "envs/smoove.yaml"
    params:
        outdir = "4_paste",
        name = PREFIX
    shell:
        """
        smoove paste --name {params.name} --outdir {params.outdir} {input}
        """