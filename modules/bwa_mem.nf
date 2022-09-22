nextflow.enable.dsl=2

process BWA_MEM {
  publishDir 'bwa_results'
  input:
  path fasta
  path fastq

  output:
  path ('*.bam'), emit: bamfile

  script:
  """
  bwa index ${fasta}
  bwa mem -t 4 ${fasta} ${fastq} > ${*.bam}
  """
}

  workflow BWA_MEM_WF {
    take:
        ch_fasta
        ch_fastq

    main:
        BWA_MEM (
          ch_fasta
          ch_fastq
          )
    emit:
        bamfile = BWA_MEM.out.bamfile
  }

  workflow {
    ch_fasta = Channel.fromPath("/mnt/c/Users/abc/Desktop/Nextflow/bwa_mem/example_human_reference.fasta")
    ch_fastq = Channel.fromPath("/mnt/c/Users/abc/Desktop/Nextflow/bwa_mem/example_human_Illumina.pe_1.fastq")
    BWA_MEM(ch_fastq)
  }
