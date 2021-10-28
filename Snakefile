configfile: "config.yaml"

from snakemake.utils import makedirs
import pandas as pd

#################################
# author: Carolina Pita Barros  #
# carolina.pitabarros@wur.nl    #
# date: August 2021             #
#################################

if "OUTDIR" in config:
    # print("\nSaving to " + config["OUTDIR"] + "\n")
    workdir: config["OUTDIR"]

makedirs("logs_slurm")

pipeline = "population-structural-var-calling-smoove" 


include: "rules/create_file_log.smk"

REFERENCE=config["REFERENCE"]
READS_DIR=config["READS_DIR"]
PREFIX = config["PREFIX"]
SAMPLES_LIST = config["SAMPLE_LIST"]
CONTIGS_IGNORE = config["CONTIGS_IGNORE"]
SPECIES = config["SPECIES"]
NUM_CHRS =  config["NUM_CHRS"]
BWA_MEM_M = config["BWA_MEM_M"]

samples_table = pd.read_csv(SAMPLES_LIST, header=None)
samples_list = list(samples_table.iloc[:,1])
SAMPLES = [os.path.splitext(x)[0] for x in samples_list]

def preprocessing(wildcards):
    '''
    If mapping has been done with bwa mem -M option, then run an extra step fo creating bam files with split and discordant reads before smoove_call
    '''
    makedirs("1_call")
    if BWA_MEM_M == "Y":
        return("1_call/{sample}.done")
    elif BWA_MEM_M == "N":
        return([])


localrules: create_file_log, simple_stats, plot_PCA

rule all:
    input:
        files_log,
        "FIGURES/" + PREFIX + ".pdf",
        "2_merged/" + PREFIX + ".smoove-counts.html",
        "6_metrics/"+ PREFIX + ".survivor.stats",
        "5_postprocessing/"+PREFIX + "_DUP_DEL_INV.vcf",
        "5_postprocessing/"+ PREFIX + "_BND.vcf",
        "5_postprocessing/"+ PREFIX + "_DUP_DEL_INV_table.tsv",


# Create comma separated list of contigs to ignore
with open(CONTIGS_IGNORE, "r") as infile:
    content =  infile.read().splitlines()
    CONTIGS =  ",".join(content)

rule split_disc_reads:
    '''
    If the mapping has been done with bwa mem -M option, using the smoove pipeline, the support for SR will always be zero. 
    This step creates the split and discordant reads beforehand, which fixes this problem
    '''
    input:
        os.path.join(READS_DIR, "{sample}.bam")
    output:
        done = touch("1_call/{sample}.done")
    message:
        'Rule {rule} processing'
    params:
        scripts_dir = os.path.join(workflow.basedir, "scripts/"),
        outdir = "1_call"
    conda:
        "envs/python2.7.yml"
    group:
        'smoove_call'
    shell:
        """
module load samtools

sname=`samtools view -H {input} | grep '^@RG' | sed "s/.*SM:\([^\\t]*\).*/\\1/g" | uniq`

python {params.scripts_dir}bamgroupreads.py -f -M -i {input} | samblaster --ignoreUnmated -M -a -e -d {params.outdir}/{wildcards.sample}.disc.sam -s {params.outdir}/{wildcards.sample}.split.sam -o /dev/null


grep -v "SAMBLASTER" {params.outdir}/{wildcards.sample}.split.sam > $sname.tmp.sam
mv $sname.tmp.sam {params.outdir}/{wildcards.sample}.split.sam
grep -v "SAMBLASTER" {params.outdir}/{wildcards.sample}.disc.sam > $sname.tmp.sam
mv $sname.tmp.sam {params.outdir}/{wildcards.sample}.disc.sam

samtools sort -@ 12 -O bam {params.outdir}/{wildcards.sample}.split.sam > {params.outdir}/$sname.split.bam
samtools sort -@ 12 -O bam {params.outdir}/{wildcards.sample}.disc.sam > {params.outdir}/$sname.disc.bam

rm {params.outdir}/{wildcards.sample}.split.sam {params.outdir}/{wildcards.sample}.disc.sam
        """


rule smoove_call:
    input:
        bam = os.path.join(READS_DIR, "{sample}.bam"),
        reference = REFERENCE,
        preprocessing = preprocessing
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
    group:
        'smoove_call'
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

smoove genotype -x \
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
        # warnings = "5_postprocessing/{prefix}.suqare.vep.vcf.gz_warnings.txt",
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
        vcf = rules.run_vep.output.vcf,
    output:
        eigenvec = "5_postprocessing/{prefix}.eigenvec",
        eigenval = "5_postprocessing/{prefix}.eigenval",
    message:
        'Rule {rule} processing'
    params:
        prefix= os.path.join("5_postprocessing",PREFIX),
        num_chrs = NUM_CHRS
    group:
        'calling'
    shell:
        """
        module load plink/1.9-180913

        plink --vcf {input.vcf} --pca --double-id --out {params.prefix} --chr-set {params.num_chrs} --allow-extra-chr --threads 8
        """

rule plot_PCA:
    input:
        eigenvec = rules.PCA.output.eigenvec,
        eigenval = rules.PCA.output.eigenval,
        sample_list = SAMPLES_LIST
    output:
        "FIGURES/{prefix}.pdf"
    message:
        'Rule {rule} processing'
    params:
        rscript = os.path.join(workflow.basedir, "scripts/basic_pca_plot.R")
    # group:
    #     'calling'
    shell:
        """
module load R/3.6.2
echo $CONDA_PREFIX
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
Rscript {params.rscript} --eigenvec={input.eigenvec} --eigenval={input.eigenval} --output={output} --sample_list={input.sample_list}
        """

rule simple_stats:
    input:
        rules.run_vep.output.vcf
    output:
        stats = "6_metrics/{prefix}.survivor.stats",
        chr_stats = "6_metrics/{prefix}.survivor.stats_CHR",
        support = "6_metrics/{prefix}.survivor.statssupport"
    message:
        'Rule {rule} processing'
    conda:
        "envs/survivor.yaml"
    params:
        tmp = "6_metrics/{prefix}.stats.temp.vcf"
    log:
        err = "logs_slurm/simple_stats_{prefix}.err",
        out = "6_metrics/{prefix}.survivor.overall.stats"
    shell:
        """
        gzip -d -c {input} > {params.tmp}
        SURVIVOR stats {params.tmp} -1 -1 -1 {output.stats} 2> {log.err} 1> {log.out}

        rm {params.tmp}
        """

rule add_depth:
    input:
        rules.run_vep.output.vcf
    output:
        dupdelinv = "5_postprocessing/{prefix}_DUP_DEL_INV.vcf",
        bnd = "5_postprocessing/{prefix}_BND.vcf"
    message:
        'Rule {rule} processing'
    params:
        script = os.path.join(workflow.basedir, "scripts/add_depth_field.py"),
        prefix = os.path.join("5_postprocessing",PREFIX)
    log:
        err = "logs_slurm/add_depth_{prefix}.err",
        out = "logs_slurm/add_depth_{prefix}.out"
    shell:
        'python {params.script} -v {input} -Q 30 -p {params.prefix} 2> {log.err} 1> {log.out}'

rule get_tsv:
    input:
        vcf = rules.add_depth.output.dupdelinv,
        samples_file = SAMPLES_LIST
    output:
        "5_postprocessing/{prefix}_DUP_DEL_INV_table.tsv"
    message:
        'Rule {rule} processing'
    params:
        script = os.path.join(workflow.basedir, "scripts/get_flat_file.py"),
    log:
        err = "logs_slurm/get_tsv_{prefix}.err",
        out = "logs_slurm/get_tsv_{prefix}.out"
    shell:
        """
        # module load python/3.9.4
        python {params.script} -v {input.vcf} -s {input.samples_file} -o {output}
        """
        