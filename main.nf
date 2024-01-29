#!/usr/bin/env nextflow
// enable dsl2
nextflow.enable.dsl=2

process DOWNOAD {
    tag { sample }
    label 'download'
    maxForks 1
    errorStrategy 'retry', 3

    publishDir "fastq/"

    input:
    tuple val(sample), val(read1), val(read2)

    output:
    tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz"), emit: fqs

    script:
    """
    wget -O ${sample}_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/$read1
    wget -O ${sample}_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/$read2
    """
}

process SNIPPY {
    tag { sample }
    //scratch 'true'

    label 'snippy'

    //publishDir "fasta/", mode: 'copy', saveAs: { filename -> "${sample}.fasta"} 

    input:
        tuple val(sample), file('reads1.fastq'), file('reads2.fastq')
    each path('ref.fa')

    output:
        tuple val(sample), file("${sample}_output/snps.consensus.subs.fa"), emit: correctedRef 
        tuple val(sample), file("${sample}_output/snps.bam"), file("${sample}_output/snps.bam.bai"), emit: bam 


    script:
    """
    snippy --cpus $task.cpus \
        --outdir ${sample}_output \
        --ref ref.fa \
        --R1 'reads1.fastq' \
        --R2 'reads2.fastq'
    """
}

process DEPTH {
    tag { sample }
    scratch 'true'
    label 'snippy'

    publishDir "pysamstats/", mode: 'copy'

    conda '/home/nick/miniconda3/envs/genericBugWorkflow'

    input:
    tuple val(sample), file("snps.bam"), file("snps.bam.bai")
    each path('ref.fa')

    output:
    tuple val(sample), path("${sample}_pysamfile.txt")

    script:
    """
    pysamstats -t variation_strand snps.bam -f ref.fa > ${sample}_pysamfile.txt
    """
}

process MASKDEPTH {
    tag { sample }
    label 'snippy'

    publishDir "masked_fasta/", pattern: '*.fasta',mode: 'copy', saveAs: { filename -> "${sample}.fasta"} 
    publishDir "masked_coverage_stats/", pattern: '*.csv', mode: 'copy', saveAs: { filename -> "${sample}.csv"} 


    input:
    tuple val(sample), path('unmasked.fa'), path('pysamfile.txt')

    output:
    tuple val(sample), path("${sample}.fasta"), emit: fasta
    tuple val(sample), path("depthStats.csv"), emit: csv


    script:
    """
    maskDepth.py -p pysamfile.txt -f unmasked.fa -o ${sample}.fasta -d 10 -mw 0.8 -r NC_011035.1:${sample}_illumina
    """
}

params.ids=''
ids = file(params.ids)
    .readLines()
    .each{it}

//ids=['ERR9633183', 'ERR9633185', 'ERR9633187', 'ERR9633188', 'ERR9633193', 'ERR9633194', 'ERR9633196', 'ERR9633198', 'ERR9633199', 'ERR9633200']
Channel.fromSRA(ids, apiKey: '9ea8db1939a1d829a43e5615d8fd4b927408')
    .map{ row -> tuple(row[0], row[1][0], row[1][1])}
    .set{illumina_fqs}

params.ref=''
Channel.fromPath(params.ref)
    .view()
    .set{reference_genome}

 workflow {

     main:

     DOWNOAD(illumina_fqs)

     SNIPPY(DOWNOAD.out.fqs, reference_genome)

     DEPTH(SNIPPY.out.bam, reference_genome)

     MASKDEPTH(SNIPPY.out.correctedRef.combine(DEPTH.out, by:0))
}

