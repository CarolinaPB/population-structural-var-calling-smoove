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
SPECIES = config["SPECIES"]

with open(SAMPLES_LIST, "r") as infile:
    content = infile.readlines()
    SAMPLES = [os.path.splitext(x)[0] for x in content]
# print(SAMPLES)

localrules: create_file_log

rule all:
    input:
        files_log,
        # expand("3_genotyped/{sample}-smoove.genotyped.vcf.gz", sample=SAMPLES),
        # expand("4_paste/{prefix}.smoove.square.vcf.gz", prefix = PREFIX),
        "FIGURES/" + PREFIX + ".pdf",
        "2_merged/" + PREFIX + ".smoove-counts.html"


test = glob_wildcards(os.path.join(READS_DIR, "/^([^.]+)/.bam"))
print(test)


with open(CONTIGS_IGNORE, "r") as infile:
    content =  infile.read().splitlines()
    CONTIGS =  ",".join(content)

rule smoove_call:
    input:
        bam = os.path.join(READS_DIR, "{sample}.bam"),
        reference = REFERENCE,
    output:
        vcf = temp("1_call/{sample}-smoove.genotyped.vcf.gz"),
        idx = temp("1_call/{sample}-smoove.genotyped.vcf.gz.csi"),
        # t1 = temp(multiext("1_call/{sample,/^([^.]+)/}.split", ".bam", ".bam.csi", ".bam.orig.bam")),
        # t2 = "1_call/{sample}.sorted-lumpy-cmd.sh",
        # t3 = temp(multiext("1_call/{sample}.disc",".bam", ".bam.csi", ".bam.orig.bam")),
        # t4 = temp("1_call/{sample}.histo")
    message:
        'Rule {rule} processing'
    params:
        outdir = "1_call",
        contigs=CONTIGS,
        scripts_dir = os.path.join(workflow.basedir, "scripts/")
    conda:
        "envs/smoove.yaml"
    shell:
        """
export PATH={params.scripts_dir}:$PATH
    
smoove call --outdir {params.outdir} \
--name {wildcards.sample} \
--fasta {input.reference} \
--duphold \
--genotype \
-p 1 \
--excludechroms {params.contigs} \
{input.bam}
        """

# TRIMMED = [os.path.splitext(f)[0] for f in SAMPLES]


rule smoove_merge:
    input:
        reference = REFERENCE,
        vcfs = expand("1_call/{sample}-smoove.genotyped.vcf.gz", sample=SAMPLES),
        idx = expand("1_call/{sample}-smoove.genotyped.vcf.gz.csi", sample=SAMPLES),
        # t1 = expand("1_call/{sample}.split{ext}", ext=[".bam", ".bam.csi", ".bam.orig.bam"], sample=SAMPLES),
        # t2 = expand("1_call/{sample}.sorted-lumpy-cmd.sh", sample=SAMPLES),
        # t3 = expand("1_call/{sample}.disc{ext}", ext=[".bam", ".bam.csi", ".bam.orig.bam"], sample=SAMPLES),
        # t4 = expand("1_call/{sample}.histo", sample=SAMPLES)
    output:
        vcf = expand("2_merged/{prefix}.sites.vcf.gz", prefix=PREFIX),
        counts = expand("2_merged/{prefix}.smoove-counts.html", prefix=PREFIX)
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
        vcf = rules.smoove_merge.output.vcf,
        bam = os.path.join(READS_DIR, "{sample}.bam"),
    output:
        vcf = "3_genotyped/{sample}-smoove.genotyped.vcf.gz",
        idx = "3_genotyped/{sample}-smoove.genotyped.vcf.gz.csi"
    message:
        'Rule {rule} processing'
    params:
        outdir = "3_genotyped",
        name = "{sample}",
        scripts_dir = os.path.join(workflow.basedir, "scripts/")
    conda:
        "envs/smoove.yaml"
    shell:
        """
export PATH={params.scripts_dir}:$PATH

smoove genotype -x -p 3 \
--name {params.name} \
--outdir {params.outdir} \
--fasta {input.reference} \
--duphold \
-p 1 \
--vcf {input.vcf} \
{input.bam}
        """

rule smoove_paste:
    input:
        expand("3_genotyped/{sample}-smoove.genotyped.vcf.gz", prefix=PREFIX, sample=SAMPLES),
    output:
        temp("4_paste/{prefix}.smoove.square.vcf.gz")
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


rule run_vep:
    input:
        rules.smoove_paste.output
    output:
        vcf = '5_postprocessing/{prefix}.smoove.square.vep.vcf.gz',
        # warnings = "results/{prefix}.vcftyper.sorted.vep.vcf.gz_warnings.txt",
        summary = "5_postprocessing/{prefix}.smoove.square.vep.vcf.gz_summary.html"
    message:
        'Rule {rule} processing'
    conda:
        "envs/vep_dependencies.yaml"
    params:
        species = SPECIES
    group:
        'calling'
    shell:
        """
module load samtools
/cm/shared/apps/SHARED/ensembl-vep/vep -i {input} \
--format vcf \
--buffer_size 5000 \
--offline \
--dir /lustre/nobackup/SHARED/cache/ \
--species {params.species} \
--vcf \
--force_overwrite \
-o {output.vcf} \
--fork 1 \
--compress_output bgzip \
--canonical \
--gene_phenotype \
--regulatory \
--numbers \
--symbol
        """

rule PCA:
    input:
        rules.run_vep.output.vcf
    output:
        eigenvec = "5_postprocessing/{prefix}.eigenvec",
        eigenval = "5_postprocessing/{prefix}.eigenval",
        pdf = "FIGURES/{prefix}.pdf"
    message:
        'Rule {rule} processing'
    params:
        prefix= os.path.join("5_postprocessing",PREFIX),
        rscript = os.path.join(workflow.basedir, "scripts/basic_pca_plot.R")
    group:
        'calling'
    shell:
        """
        module load R/3.6.2
        module load plink/1.9-180913

        plink --vcf {input} --pca --double-id --out {params.prefix} --chr-set 38 --allow-extra-chr --threads 8
        Rscript {params.rscript} --eigenvec={output.eigenvec} --eigenval={output.eigenval} --output={output.pdf}
        """