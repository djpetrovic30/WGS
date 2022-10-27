nextflow.enable.dsl=2

process snpEff {
  publishDir params.annvcf

  input:
  path vcf

  output:
  path "ann20.vcf", emit: annvcf

  script:
  """
  java -Xmx8g -jar /home/biodocker/bin/snpEff/snpEff.jar ann -v GRCh37.75 ${vcf} > ann20.vcf
  """
  }

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow snpEff_wf{
  take:
  ch_vcf
  main:
  snpEff(vcf)
  emit:
  annvcf
}

workflow {
  ch_vcf = Channel.fromPath(params.vcf)
  snpEff(ch_vcf)
}
