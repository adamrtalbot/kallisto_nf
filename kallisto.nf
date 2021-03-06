#!/usr/bin/env nextflow

/// Input variables

params.reference = "*.fa"
params.kmer = 31
params.reads = "*{1,2}.fastq.gz"
params.boostraps = 100
params.output = "output"
params.index_memory = "4G"
params.quant_memory = "4G"
params.tximport_memory = "4G"
params.threads = 4
params.time = "1h"

/// Reference to channel for making reference
/// Set read pairs as channel, based on pattern match specified in input.

cdna_fasta = file(params.reference)

Channel
  .fromFilePairs( params.reads )
  .ifEmpty { error "Cannot find reads matching: ${params.reads}" }
  .set { read_pairs }


/// Make reference from cDNA sequences

process makeReference {

  tag "Making cDNA reference index"

  module 'kallisto'

  time "${params.time}"
  clusterOptions "-cwd -V -S /bin/bash -l h_vmem=${params.index_memory}"

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

  publishDir = [path: params.output, mode: 'copy']

  tag "read: ${name}"

  module 'kallisto'

  time "${params.time}"
  clusterOptions "-cwd -V -S /bin/bash -pe smp ${params.threads} -l h_vmem=${params.quant_memory}"

  input:
  file idx from kallisto_index
  set name, file(fastq_pair) from read_pairs

  output:
  file "${name}" into kallist_out

  """
  kallisto quant -i ${idx} -b ${params.boostraps} -t ${params.threads} -o ${name} ${fastq_pair}
  """

}


/// Collate all kallisto .hd5 files and arrange into gene count table using tximport

process  tximport {

  publishDir = [path: params.output, mode: 'copy']

  tag "Collecting count data"

  module 'R'

  time "${params.time}"
  clusterOptions "-cwd -V -S /bin/bash -l h_vmem=${params.tximport_memory}"

  input:
  file kallisto_abundance from kallist_out.collect()

  output:
  file "count_data.csv" into kallisto_gathered

  """
  tximport_kallisto.R --vanilla "${workflow.projectDir}/${params.output}"
  """

}
