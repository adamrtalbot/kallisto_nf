#!/usr/bin/env nextflow

/// Input variables

params.reference = "*.fa"
params.kmer = 31
params.reads = "*{1,2}.fq.gz"
params.boostraps = 100

/// Reference to channel for making reference
/// Set read pairs as channel, based on pattern match specified in input.

cdna_fasta = Channel.fromPath(params.reference)
read_pairs = Channel.fromFilePairs(params.reads)

/// Make reference from cDNA sequences

process makeReference {

  input:
  file 'fasta' from cdna_fasta

  output:
  file 'kallisto.idx' into kallisto_index

  """
  kallisto index -i kallisto.idx -k ${params.kmer} fasta
  """

}

/// Run Kallisto Quant on all read file pairs

process quantReads {

  input:
  file index from kallisto_index
  file fastq_pair from read_pairs

  output:
  file mapping into kallist_out_dirs

  """
  kallisto quant -i index -b ${params.boostraps} -t 1 -o mapping fastq_pair
  """

}
