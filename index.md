---
layout: page
title: Population SV calling with smoove
group: navigation
---

[Link to the repository](https://github.com/CarolinaPB/population-structural-var-calling-smoove/tree/single_run)

## First follow the instructions here:
[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT
This is a pipeline to perform structural variant calling in a population using Smoove. It also runs VEP and performs PCA. 
In addition to the VCF with the SVs, you also get a .tsv file with some summarized information on the SVs: it included allele frequency per population, as well as VEP annotation and depth fold change as described in [duphold](https://github.com/brentp/duphold:
> DHBFC: fold-change for the variant depth relative to bins in the genome with similar GC-content.  
> DHFFC: fold-change for the variant depth relative to Flanking regions.


#### Tools used:
- Smoove - SV calling
- VEP - determines the effect of the variants
- Plink - perform PCA
- R - plot PCA
- SURVIVOR - basic SV stats
- Python - add depth to vcf and create final table 
  - PyVcf (you might need to install it)


| ![DAG](https://github.com/CarolinaPB/population-structural-var-calling-smoove/blob/single_run/dag.png) |
|:--:|
|*Pipeline workflow* |


### Edit config.yaml with the paths to your files
```
OUTDIR: /path/to/output 
READS_DIR: /path/to/reads/ # don't add the reads files, just the directory where they are
SAMPLE_LIST: /path/to/file
REFERENCE: /path/to/assembly
CONTIGS_IGNORE: /path/to/file
SPECIES: <species_name>
PREFIX: <output name>
```

- OUTDIR - directory where snakemake will run and where the results will be written to
- READS_DIR - path to the directory that contains the reads
- SAMPLE_LIST - three column csv with the sample name, name of the bam files to use in the second column and the name of the corresponding population on the third column. These bams should all be in the same directory (READS_DIR)
- Example: 
> sample1,sample1.bam,Pop1   
> sample2,sample2.bam,Pop1   
> sample3,sample3.bam,Pop2   
> sample4,sample4.bam,Pop2  
- REFERENCE - path to the assembly file
- CONTIGS_IGNORE - contigs to be excluded from SV calling (usually the small contigs)
- SPECIES - species name to be used for VEP
- PREFIX - prefix for the created files

If you want the results to be written to this directory (not to a new directory), comment out or remove
```
OUTDIR: /path/to/outdir
```

## ADDITIONAL SET UP
This pipeline uses VEP in offline mode, which increases performance. In order to use it in this mode, the cache for the species used needs to be installed:
#### For people using WUR's Anunna:
Check if the cache file for your species already exist in `/lustre/nobackup/SHARED/cache/`. If it doesn't, create it with

```
/usr/bin/perl /cm/shared/apps/SHARED/ensembl-vep/INSTALL.pl --CACHEDIR /lustre/nobackup/SHARED/cache/ --AUTO c -n --SPECIES <species>
```
When multiple assemblies are found you need to run it again with `--ASSEMBLY <assembly name>`, where "assembly name" is the name of the assembly you want to use.

#### For those not from WUR:
You can install VEP with 
```
conda install -c bioconda ensembl-vep
```
and install the cache with 
```
vep_install --CACHEDIR <where/to/install/cache> --AUTO c -n --SPECIES <species>
```
When multiple assemblies are found you need to run it again with `--ASSEMBLY <assembly name>`, where "assembly name" is the name of the assembly you want to use.

In the Snakefile, in rule `run_vep`, replace `/cm/shared/apps/SHARED/ensembl-vep/vep` with `vep`

## RESULTS
* **<run_date>_files.txt** Dated file with an overview of the files used to run the pipeline (for documentation purposes)
* **2_merged** 
  * {prefix}.smoove-counts.html - shows a summary of the number of reads before and after filtering 
* **5_postprocessing** directory that contains the final VCF file containing the structural variants found. This file has been annotated with VEP
  * {prefix}.smoove.square.vep.vcf.gz - Final VCF
  * {prefix}.smoove.square.vep.vcf.gz_summary - statistics from VEP
  * {prefix}.nosex, {prefix}.log, {prefix}.eigenvec, {prefix}.eigenval - output files from the PCA
  * {prefix}_DUP_DEL_INV_table.tsv - table with the most important information extracted from the VCF. Contains information about the SV, allele frequency for each population, VEP annotation and depth information
* **6_metrics** directory that contains general stats about the number of SVs found
