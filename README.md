# kallisto_nf
Nextflow pipeline for Kallisto pseudoalignment of RNA-Seq data


## Description

This is a nextflow pipeline for running Kallisto against RNA-Seq reads. It will:
1. Generate a kallisto index from a fasta file of cDNA sequences (specified by --reference)
2. Run pseudoaligment and quantification using kallisto quant (specified in --reads)
3. Import and aggregate all pseudocounts with tximport in R (using the Rscript file in bin/tximport_kallisto.R)

## Usage

To run on a server, you need enter this as a command:

```java
nextflow run adamrtalbot/kallisto_nf \
--reference path/to/reference/cDNA.fasta \
--reads 'path/to/reads/reads_{1,2}.fastq.gz'
-profile sge \
-bg
```

To break this down:

```java
nextflow run adamrtalbot/kallisto_nf
```

Runs the pipeline from this repository.

```java
--reference path/to/reference/cDNA.fasta
```

Path to reference, composed of a fasta file of all cDNA (transcript) sequences.

```java
--reads 'path/to/reads/reads_{1,2}.fastq.gz'
```

Specify path to reads. Use the two variables in the curly brackets (e.g. ```{1,2}```) to specify forward and reverse read. Note it needs to be in quotation marks, this is for Nextflow to parse the forward and reverse reads properly.

```java
-profile sge
```

Specifies this is to be run as if on a cluster, using the Sun Grid Engine. This can be specified to ```-profile standard``` to run on your local machine. Change the nextflow.config file to add more profiles that fit your specific use case.

```java
-bg
```

Run in background. You may have to hit enter to exit the nextflow pipeline.

Add ```-resume``` to resume a previously ran (and presumably failed) analysis.

## Optional

There are a few optional parameters you can change, they have default options specified below:

```java
// Number of cores
--threads // = 4
// Number of Kallisto quantification bootstraps
--boostraps // = 100
// Output folder
--output // = output
// Amount of RAM supplied to each process
--memory // = 4G (4gb memory)
// Kmer size for Kallisto index
--kmers // = 31

Full configuration format is standard Nextflow, as such you may need to add more profiles etc. to work with your particular HPC set up.

## Logs and reports etc

I have no added extra logging or any tools. If you would like more, use the built in tools available here: https://www.nextflow.io/docs/latest/tracing.html
