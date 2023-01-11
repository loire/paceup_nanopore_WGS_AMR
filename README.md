# paceup_nanopore_WGS_AMR
snakemake pipeline dedicated to analysis of WGS analysis of resistant bacteria 

## Install

## steps 

from raw fastq:

* cleaning with FiltLong
* assembly with Flye
* Busco analysis of assemblies
* Polishing with medaka
* AMR gene search with Abricate
* Mobile element search with Mob suite
* Annotation with prokka
* core genome alignment with roary
* Phylogenomic tree with FastTree

## Usage
`snakemake -s Nanopore_bacterial_WGS_Assembly.snakefile -r --cores 1`

## Credits


