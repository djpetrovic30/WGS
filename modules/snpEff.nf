nextflow.enable.dsl=2

process snpEff {
  publishDir params.annvcf

  input:
  path vcf
  path config

  output:
  path "ann20.vcf", emit: annvcf

  script:
  """
  java -Xmx8g -jar /home/biodocker/bin/snpEff/snpEff.jar -c ${config} ann -v humanchr20 ${vcf} > ann20.vcf
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
  vcf
  config
  main:
  snpEff(vcf, config)
  emit:
  annvcf
}

workflow {
  ch_vcf = Channel.fromPath(params.vcf)
  ch_config = Channel.fromPath(params.config)
  snpEff(ch_vcf, ch_config)
}
