nextflow.enable.dsl=2

process Samtools_sort {
  input:
  path(bamfiles)

  output:
  path "bamsorted", emit: bamsorted

  script:
  """
  samtools sort ${bamfiles} > bamsorted
  """
}

process Qualimap_bamqc {
  // tag "${sample_id}"
  publishDir params.qualimap_bamqc
  input:
  path(bamsorted)

  // tuple val(sample_id), path(bamfiles)

  // output:
  // path("*"), emit: bamqc
  // tuple val(sample_id), path("${sample_id}_bamqc.bam"), emit: bamqc_bam

  script:
  """
  qualimap bamqc -bam ${bamsorted} -outdir $params.qualimap_bamqc
  """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
workflow Samtools_sort_wf {
  take:
    bamfiles
  main:
    Samtools_sort(bamfiles)
  emit:
    bamsorted = Samtools_sort.out.bamsorted
}
workflow Qualimap_bamqc_wf {
  take:
    bamsorted
  main:
    Qualimap_bamqc(bamsorted)
  emit:
    bamqc = Qualimap_bamqc.out.bamqc
  }

  workflow {
    ch_bamfiles = Channel.fromPath(params.bamfiles)
    Samtools_sort(ch_bamfiles)
    Qualimap_bamqc(Samtools_sort.out)
  }
