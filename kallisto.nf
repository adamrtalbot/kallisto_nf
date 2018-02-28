#!/usr/bin/env nextflow

/// Input variables

params.reference = "*.fa"
params.kmer = 31
params.reads = "*{1,2}.fastq.gz"
params.boostraps = 100

/// Reference to channel for making reference
/// Set read pairs as channel, based on pattern match specified in input.

cdna_fasta = file(params.reference)

Channel
  .fromFilePairs( params.reads, size: -1 )
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

  input:
  file idx from kallisto_index
  set val(name), file(fastq_pair) from read_pairs

  output:
  file "${name}_out" into kallist_out

  """
  kallisto quant -i ${idx} -b ${params.boostraps} -o ${name}_out ${fastq_pair}
  """

}
