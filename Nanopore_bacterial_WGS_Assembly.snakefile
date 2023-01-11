import glob
import re
import sys
from os.path import join 

###############################################################################
# Pipeline to analyze minion data from resistant bacteria isolates
############################################################################## 

#### Please edit your config file (Nanopore_bacterial_WGS_Assembly_params.yml) and keep it in the same directory as this script when you execute it

#### Commande exemple:
### TODO

configfile: "Nanopore_bacterial_WGS_Assembly_params.yml"

datadir=config["datadir"]
ext=config["ext"]
outputdir=config["outputdir"]
minreadslength=config["minReadsLength"]
bestqualreadspercentkeep=config["BestQualReadsPercentKeep"]
threadsflye=config["threadsFlye"]


################################################################################

##### Get sample lists with wildcards

SAMPLES, = glob_wildcards(datadir+"/{sample}"+ext)

##### Define helper functions for standard output during execution
def message(mes):
	sys.stderr.write("|---- " + mes + "\n")

def errormes(mes):
	sys.stderr.write("| ERROR ----" + mes + "\n")

####### Print info about samples:

NBSAMPLES = len(SAMPLES)
message(str(NBSAMPLES)+" samples  will be analysed:")
for i in SAMPLES:
	message(str(i))

 
##### Define Rules !

rule All:
	input:
		expand(outputdir+"/{smp}/FilteredReads/Filt_{smp}.fastq",smp=SAMPLES),
		expand(outputdir+"/{smp}/FlyeAssemblies/Flye_{smp}.fasta",smp=SAMPLES),
		expand(outputdir+"/{smp}/PolishedAssemblies/{smp}.fasta",smp=SAMPLES),
		expand(outputdir+"/{smp}/ResistanceAbricate/{smp}_abricate.txt",smp=SAMPLES),
		expand(outputdir+"/{smp}/BuscoOutput/{smp}_busco.txt",smp=SAMPLES),
		expand(outputdir+"/{smp}/BuscoOutput/{smp}_busco_parsed.txt",smp=SAMPLES),
		expand(outputdir+"/{smp}/Mob_recon/contig_report.txt",smp=SAMPLES),
		expand(outputdir+"/{smp}/Prokka/{smp}.gff",smp=SAMPLES),
		outputdir+"/roary/core_gene_alignment.aln",
		outputdir+"/FastTree/tree.newick"

rule Filter:
	input:
		datadir+"/{smp}"+ext
	output:
		outputdir+"/{smp}/FilteredReads/Filt_{smp}.fastq"
	params:
		minl = config["minReadsLength"],
		bqual = config["BestQualReadsPercentKeep"]
	shell:
		"""
		filtlong --min_length {params.minl} --keep_percent {params.bqual} {input} > {output}
		"""
rule Flye:
	input:
		outputdir+"/{smp}/FilteredReads/Filt_{smp}.fastq" 
	output:
		outputdir+"/{smp}/FlyeAssemblies/Flye_{smp}.fasta"
	params:
		outfly = outputdir+"/{smp}/FlyeAssemblies/",
		thread = config["threadsFlye"]
	shadow: "shallow"
	shell:
		"""
		mkdir -p {params.outfly}
		flye --nano-hq {input} --threads {params.thread} --out-dir {params.outfly}  && cp {params.outfly}/assembly.fasta {output}
		"""
rule Polish:
	input:
		reads = outputdir+"/{smp}/FilteredReads/Filt_{smp}.fastq",
		assembly = outputdir+"/{smp}/FlyeAssemblies/Flye_{smp}.fasta"
	output:
		outputdir+"/{smp}/PolishedAssemblies/{smp}.fasta"
	shadow: "shallow"
	params:
		medaCPU = config["threadsMedaka"],
		medaMod = config["medakamodel"],
		medaBatch = config["batchMedaka"],
		outdir = outputdir+"/{smp}/PolishedAssemblies"
	shell:
		"""
		medaka_consensus -i {input.reads}  -d {input.assembly} -t {params.medaCPU} -f -m {params.medaMod} -b {params.medaBatch}
		cp -r medaka {params.outdir}
		cp medaka/consensus.fasta {output} 
		"""		


rule abricate:
	input:
		outputdir+"/{smp}/PolishedAssemblies/{smp}.fasta"
	output:
		outputdir+"/{smp}/ResistanceAbricate/{smp}_abricate.txt"
	shell:
		"""
		abricate {input} > {output}
		"""
rule busco:
	input:
		outputdir+"/{smp}/PolishedAssemblies/{smp}.fasta"
	output:
		outputdir+"/{smp}/BuscoOutput/{smp}_busco.txt"
	shadow: 
		"shallow"
	params:
		buscoL = config["buscoLineage"]
	shell:
		"""
		busco -m genome -i {input} -o tmp -l {params.buscoL}
		cp tmp/short_summary.specific*.txt {output}
		"""

rule parsebusco:
	input: 
		outputdir+"/{smp}/BuscoOutput/{smp}_busco.txt"
	output:
		outputdir+"/{smp}/BuscoOutput/{smp}_busco_parsed.txt"
	shell:
		"""
		awk -F"\\t" 'NR>9 && NR<16{{head=head$3"\\t" ; val=val$2"\\t"}}END{{print head;print val}}' {input} > {output}
		"""

rule mob_recon:
	input:
		outputdir+"/{smp}/PolishedAssemblies/{smp}.fasta"
	output:
		outputdir+"/{smp}/Mob_recon/contig_report.txt"
	shadow:
		"shallow"
	params:
		outdir=outputdir+"/{smp}/Mob_recon/",
		mobthreads=config["threadsMob"]
	shell:
		"""
		mob_recon --force -i {input} -o {params.outdir} -n {params.mobthreads}
		"""
rule prokka:
	input:
		outputdir+"/{smp}/PolishedAssemblies/{smp}.fasta"
	output:
		outputdir+"/{smp}/Prokka/{smp}.gff"
	params:
		outdir=outputdir+"/{smp}/Prokka",
		prefix="{smp}"
	shell:
		"""
		prokka --outdir {params.outdir} --prefix {params.prefix} --force {input}
		"""
rule roary:
	input: 
		expand(outputdir+"/{smp}/Prokka/{smp}.gff",smp=SAMPLES)		
	output:
		outputdir+"/roary/core_gene_alignment.aln"
	params:
		outdir=outputdir+"/roary",
		threads=config["threadsRoary"]
		
	shell:
		"""
		roary -p {params.threads} -e --mafft -f {params.outdir} {input}
		"""

rule fasttree:
	input:
		outputdir+"/roary/core_gene_alignment.aln"
	output:
		outputdir+"/FastTree/tree.newick"
	shell:
		"""
		FastTree -gtr -nt {input} > {output}
		"""


#
#rule MakeTableGraphs:
#	input:
#		expand() ...
#
#	shell:
#		"""
#		Rscript get_summary_stat.R {input}
#		touch {output}
#		"""
#

