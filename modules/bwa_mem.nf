nextflow.enable.dsl=2

process BWA_MEM {
  publishDir 'bwa_results'
  input:
  path fasta
  path fastq

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

  workflow BWA_MEM_WF {
    take:
        ch_fasta
        ch_fastq

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
    emit:
        bamfile = BWA_MEM.out.bamfile

  }

  workflow {
    ch_fasta = Channel.fromPath(params.fasta)
    ch_fastq = Channel.fromPath(params.fastq)
     BWA_MEM_WF(ch_fastq, ch_fasta)
  }
