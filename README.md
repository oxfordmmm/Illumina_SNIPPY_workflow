# Illumina Workflow README

This repository contains the `main.nf` file which is used as the workflow for the Nextflow pipeline. 

## Running the Workflow

To run the workflow, use the following command:

```bash
nextflow run ~/soft/illumina_workflow/main.nf \
        --ref NZ_CP007601.1.fasta \
        --ids list_of_accessions.txt \
        -resume
```
## Parameters
*`--ref`: This is the reference genome file in FASTA format. In the example above, `NZ_CP007601.1.fasta` is used.
*`--ids`: This is a text file containing the list of accession numbers. In the example above, `list_of_accessions.txt` is used.
*`-resume`: This flag allows Nextflow to resume the pipeline from where it last left off.

## Dependencies
Dependencies are managed by the nextflow.config file, which points to the Conda environment YAML files located in the `env/` directory. This ensures that all required software and libraries are installed and available for the workflow.


