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
  input:
  path bamfile

  output:
  tuple val, path('example_sorted.bam'), path('example_sorted.bai'), emit: bamsorted

  script:
  """
  samtools sort ${bamfile} example_sorted.bam
  samtools index example_sorted.bam > example_sorted.bai
  """
}

process MarkDuplicates {
  input:
  path 'example_sorted.bam'

  output:
  path 'marked_duplicates.bam'

  script:
  """
  java -jar /root/gatk-package-4.2.6.1/gatk-package-4.2.6.1-local.jar MarkDuplicates \
  -I 'exanple.bam' \
  -O 'marked_duplicates.bam' \
  -M 'marked_dup_metrics.txt'
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
    emit:
        bamfile = BWA_MEM.out.bamfile
  }

  workflow {
    ch_fasta = Channel.fromPath("/mnt/c/Users/abc/Desktop/Nextflow/bwa_mem/example_human_reference.fasta")
    ch_fastq = Channel.fromPath("/mnt/c/Users/abc/Desktop/Nextflow/bwa_mem/example_human_Illumina.pe_1.fastq")
     BWA_MEM(ch_fastq, ch_fasta)
  }
