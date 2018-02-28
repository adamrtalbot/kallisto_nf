#!/usr/bin/env nextflow

/// Input variables

params.reference = "*.fa"
params.kmer = 31
params.reads = "*{1,2}.fastq.gz"
params.boostraps = 100
params.output = "output"

/// Reference to channel for making reference
/// Set read pairs as channel, based on pattern match specified in input.

cdna_fasta = file(params.reference)

Channel
  .fromFilePairs( params.reads )
  .ifEmpty { error "Cannot find reads matching: ${params.reads}" }
  .set { read_pairs }


/// Make reference from cDNA sequences

process makeReference {

  input:
  file cdna_fasta

  output:
  file "kallisto.idx" into kallisto_index

  script:
  """
  kallisto index -i kallisto.idx -k ${params.kmer} ${cdna_fasta}
  """

}

/// Run Kallisto Quant on all read file pairs

process quantReads {

  publishDir "$output", mode: 'copy'

  tag "read: ${name}"

  input:
  file idx from kallisto_index
  set name, file(fastq_pair) from read_pairs

  output:
  file "${name}" into kallist_out

  """
  kallisto quant -i ${idx} -b ${params.boostraps} -o ${name} ${fastq_pair}
  """

}


/// Collate all kallisto .hd5 files and arrange into gene count table using tximport
