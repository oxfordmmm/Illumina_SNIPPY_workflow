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

process REFMASK {
    label 'blast'
    publishDir "refmask", mode: "copy"

    input:
	path('refFasta')

	output:
    path("rpt_mask.gz"), emit: masked_ref
    

	script:
	"""
    makeblastdb -dbtype nucl -in refFasta
    genRefMask.py -r refFasta -m 200 -p 95
    bgzip -c refFasta.rpt.regions > rpt_mask.gz
	echo '##INFO=<ID=RPT,Number=1,Type=Integer,Description="Flag for variant in repetitive region">' > rpt_mask.hdr
	tabix -s1 -b2 -e3 rpt_mask.gz
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
    tuple val(sample), path('unmasked.fa'), path('pysamfile.txt'), path('masked_ref.bed.gz')

    output:
    path("${sample}.fasta"), emit: fasta
    tuple val(sample), path("depthStats.csv"), emit: csv


    script:
    """
    maskDepth.py -p pysamfile.txt \
        -f unmasked.fa \
        -o ${sample}.fasta \
        -d 10 -mw 0.8 \
        -b masked_ref.bed.gz \
        -r NZ_CP007601.1:${sample}_illumina
    """
}

process SNP_DISTS {
    tag { "create a pairwise SNP matrix" }

    label 'snpdist'
    publishDir "snp-dists", mode: "copy"

    input:
    path(nonrec)

    output:
    path("*")
    //path("${nonrec}.snp-dists.tsv"), emit: snp_matrix_tsv
    //path("${nonrec}.snp-dists.csv"), emit: snp_matrix_csv
    //path("${nonrec}.snp-dists_molten.csv"), emit: snp_matrix_molten

    script:
    """
    snp-dists ${nonrec} > ${nonrec}.snp-dists.tsv
    snp-dists -c ${nonrec} > ${nonrec}.snp-dists.csv
    snp-dists -c ${nonrec} -m > ${nonrec}.snp-dists_molten.csv
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

     REFMASK(reference_genome)

     SNIPPY(DOWNOAD.out.fqs, reference_genome)

     DEPTH(SNIPPY.out.bam, reference_genome)

     MASKDEPTH(SNIPPY.out.correctedRef.combine(DEPTH.out, by:0).combine(REFMASK.out.masked_ref))

    //comine masked fasta output into a single file
     MASKDEPTH.out.fasta.view()
        .collectFile(name: 'all_masked.fasta')
        .set{masked}  

    SNP_DISTS(masked)

 }

