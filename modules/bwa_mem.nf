nextflow.enable.dsl=2

process BWA_MEM {
  tag "${sample_id}"
  publishDir 'bwa_results'
  input:
  path fasta
  tuple val(sample_id), path(fastq)

  output:
  tuple val(sample_id), path("${sample_id}.bam"), emit: bamfile

  script:
  """
  bwa index ${fasta}
  bwa mem -t 4 ${fasta} ${fastq} > ${sample_id}.bam
  """
}

process Samtools_sort {
  publishDir 'bwa_results'
  input:
  tuple val(sample_id), path(bamfile)

  output:
  tuple val(sample_id), path("${sample_id}.sorted_bam"), emit: bamsorted
  tuple val(sample_id), path("${sample_id}.sorted_bai"), emit: baisorted

  script:
  """
  samtools sort ${bamfile} > ${sample_id}.sorted_bam
  samtools index ${sample_id}.sorted_bam > ${sample_id}.sorted_bai
  """
}

process MarkDuplicates {
  publishDir 'bwa_results'
  input:
  tuple val(sample_id), path(bamsorted)
  tuple val(sample_id), path(baisorted)

  output:
  tuple val(sample_id), path("${sample_id}_duplicates.bam"), emit: deduplicate_bam
  path "${sample_id}_metrics.txt"

  script:
  """
  gatk AddOrReplaceReadGroups I=${bamsorted} O=output.bam RGLB=lib1 RGPL=ILLUMINA RGPU=50 RGSM=TAAGGCGA
  gatk MarkDuplicates \
  -I output.bam \
  -O ${sample_id}_duplicates.bam \
  -M ${sample_id}_metrics.txt
  """
}

process BaseRecalibrator {
  publishDir 'bwa_results'
  input:
  tuple val(sample_id), path(deduplicate_bam)
  path fasta
  path dbsnp
  path dbsnpindex
  path index
  path fasta_dict
  output:
  path("recaldata.table"), emit: recalibration_table

  script:
  """
  gatk BaseRecalibrator \
  -I ${deduplicate_bam} \
  -R ${fasta} \
  --known-sites ${dbsnp} \
  --output recaldata.table
  """
}

process ApplyBQSR {
  publishDir 'bwa_results'
  input:
  path fasta
  path index
  path fasta_dict
  tuple val(sample_id), path(deduplicate_bam)
  path recalibration_table

  output:
  path "${sample_id}_recal.bam", emit: recalbam

  script:
  """
  gatk ApplyBQSR \
  -R ${fasta} \
  -I ${deduplicate_bam} \
  --bqsr-recal-file ${recalibration_table} \
  -O ${sample_id}_recal.bam
  """

}

  workflow BWA_MEM_wf {
    take:
        fasta
        fastq

    main:
        BWA_MEM (fasta, fastq)
    emit:
        bamfile = BWA_MEM.out.bamfile
        }

  workflow Samtools_sort_wf {

    take:
        bamfile
    main:
        Samtools_sort(bamfile)
    emit:
        bamsorted = Samtools_sort.out.bamsorted
        baisorted = Samtools_sort.out.baisorted

    }

  workflow MarkDuplicates_wf {
    take:
      bamsorted
      baisorted
    main:
      MarkDuplicates(bamsorted, baisorted)
    emit:
      deduplicate_bam = MarkDuplicates.out.deduplicate_bam
  }

  workflow BaseRecalibrator_wf {
    take:
      deduplicate_bam
      fasta
      dbsnp
      dbsnpindex
      index
      fasta_dict
    main:
      BaseRecalibrator(deduplicate_bam, fasta, dbsnp, dbsnpindex, index, fasta_dict)
    emit:
      recalibration_table = BaseRecalibrator.out.recalibration_table
  }

  workflow ApplyBQSR_wf {
    take:
      fasta
      index
      fasta_dict
      deduplicate_bam
      recalibration_table

    main:
      ApplyBQSR(fasta, index, fasta_dict, deduplicate_bam, recalibration_table)

    emit:
      recalbam = ApplyBQSR.out.recalbam

  }

  workflow {
    ch_fasta = Channel.fromPath(params.fasta)
    ch_fastq = Channel.fromFilePairs(params.fastq+"/*_{1,2}.fastq")
    ch_dbsnp = Channel.fromPath(params.dbsnp)
    ch_dbsnpindex = Channel.fromPath(params.dbsnpindex)
    ch_index = Channel.fromPath(params.index)
    ch_fasta_dict = Channel.fromPath(params.fasta_dict)
    BWA_MEM_wf(ch_fasta, ch_fastq)
    Samtools_sort_wf(BWA_MEM_wf.out)
    MarkDuplicates_wf(Samtools_sort_wf.out)
    BaseRecalibrator_wf(MarkDuplicates_wf.out, ch_fasta, ch_dbsnp, ch_dbsnpindex, ch_index, ch_fasta_dict)
    ApplyBQSR_wf(ch_fasta, ch_index, ch_fasta_dict, MarkDuplicates_wf.out, BaseRecalibrator_wf.out)
  }
