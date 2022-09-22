#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SAMTOOLS_VIEW {
 publishDir 'results'
  input:
  path bamfile

  output:
  path ("*.sam"), emit: samfile

 script:

 """
  samtools view -h example.bam > example.bam.sam
 """
}


workflow SAMTOOLS_VIEW_WF {
  take:
      ch_bamfile
  main:
      SAMTOOLS_VIEW (
         ch_bamfile
                  )

 emit:
     samfile = SAMTOOLS_VIEW.out.samfile
}

workflow {
  ch_bamfile = Channel.fromPath("/mnt/c/Users/abc/Desktop/Nextflow/samtools/docker/example.bam")
  SAMTOOLS_VIEW(ch_bamfile)
}
