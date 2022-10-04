nextflow.enable.dsl=2

process BWA_MEM {
//  tag "${sample_id}"
  publishDir 'bwa_results'
  input:
  path fasta
  path fastq
  // path fastq

  output:
  path ('example.bam'), emit: bamfile

  script:
  """
  bwa index ${fasta}
  bwa mem -t 4 ${fasta} ${fastq} > example.bam
  """
}

process Samtools_sort {
  publishDir 'bwa_results'
  input:
  path bamfile

  output:
  path('example_sorted.bam'), emit: bamsorted
  path('example_sorted.bai'), emit: baisorted

  script:
  """
  samtools sort ${bamfile} > example_sorted.bam
  samtools index example_sorted.bam > example_sorted.bai
  """
}

process MarkDuplicates {
  publishDir 'bwa_results'
  input:
  path bamsorted
  path baisorted

  output:
  path 'marked_duplicates.bam', emit: deduplicate_bam
  path 'marked_dup_metrics.txt', emit: deduplicate_metrics

  script:
  """
  gatk MarkDuplicates \
  -I ${bamsorted} \
  -O marked_duplicates.bam \
  -M marked_dup_metrics.txt
  """
}

process BaseRecalibrator {
  publishDir 'bwa_results'
  input:
  path deduplicate_bam
  path fasta
  path index
  path dbsnp
  path fasta_dict
  output:
  path 'recal_data.table', emit: recalibration_table

  script:
  """
  gatk BaseRecalibrator \
  -I ${deduplicate_bam} \
  -R ${fasta} \
  --known-sites ${dbsnp} \
  --output recal_data.table
  """
}

// process ApplyBQSR {
//   publishDir 'bwa_results'
//   input:
//   path recalibration.table
//
//   output:
//   path 'recal.bam', emit: recalbam
//
//   script:
//   """
//   gatk ApplyBQSR \
//   -R ${fasta} \
//   -I ${deduplicate_bam} \
//   --bqsr-recal-file ${recalibration.table} \
//   -O recal.bam
//   """
//
// }

  workflow BWA_MEM_WF {
    take:
        ch_fasta
        ch_fastq
        ch_index
        ch_dbsnp
        ch_fasta_dict

    main:
        BWA_MEM (
          ch_fasta,
          ch_fastq
          )

      ch_bam = BWA_MEM.out.bamfile
          Samtools_sort (
            ch_bam
            )
      ch_bamsorted = Samtools_sort.out.bamsorted
      ch_baisorted = Samtools_sort.out.baisorted
        MarkDuplicates (
           ch_bamsorted, ch_baisorted
          )

      ch_deduplicate_bam = MarkDuplicates.out.deduplicate_bam
         BaseRecalibrator (
           ch_deduplicate_bam, ch_fasta, ch_index, ch_dbsnp, ch_fasta_dict
           )
    emit:
        bamfile = BWA_MEM.out.bamfile

  }

  workflow {
    ch_fasta = Channel.fromPath(params.fasta)
    ch_fastq = Channel.fromPath(params.fastq)
    ch_dbsnp = Channel.fromPath(params.dbsnp)
    ch_index = Channel.fromPath(params.index)
    ch_fasta_dict = Channel.fromPath(params.fasta_dict)
    BWA_MEM_WF(ch_fasta, ch_fastq, ch_index, ch_dbsnp, ch_fasta_dict)
  }
