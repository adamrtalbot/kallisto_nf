#!/usr/bin/env nextflow

/// Input variables

params.reference = "*.fa"
params.kmer = 31

/// Add input variable to channels

cdna_fasta = Channel.fromPath(params.reference)


/// Make reference from cDNA sequences

process makeReference {

  input:
  file 'fasta' from cdna_fasta

  output:
  file 'kallisto.idx' into kallisto_index

  """
  cat $fasta > transcriptome
  kallisto index -i kallisto_index -k ${params.kmer} transcriptome
  """

}
