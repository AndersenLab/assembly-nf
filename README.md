Nextflow pipeline for PacBio sequencing quality control, *de novo* genome assembly, and de-contamination of non-*Nematoda* DNA.

<img alt="assembly-nf/workflow metro map" src="https://github.com/AndersenLab/assembly-nf/blob/main/workflow_animation.svg">

Workflow diagram generated using [nf-metro](https://github.com/seqeralabs/nf-metro).

By default, the pipeline will assemble genomes in "inbred" mode (i.e., will not purge duplicate haplotigs) for hifiasm.

This pipeline can produce assemblies using just HiFi seq data (e.g., ```main_wBlobToolKit.nf```) or with HiFi and ONT data (e.g., ```main_wBlobToolKit_andONT.nf```). If ```main_wBlob_andONT.nf``` is used, the user must add a column to the supplied ```--sample_sheet``` with header ```ont_path``` and the absolute path pointing to the ONT seq data.

This pipeline expects unaligned HiFi BAM files and gzipped ONT FASTQ files. We recommend filtering ONT reads for >= 10 kb in length and Phred >= 15 prior to supplying to ```assembly-nf```.

## --source
If running with --source umd, then provide the folder in ```/vast/eande106/data/transfer/raw``` where sequencing data has been deposited (e.g. ```--raw_dir 20250314_PacBio```). When running in ```--source umd```, the output directory will be set to ```"${raw_dir}-assembly"``` even if you specify an ```--outdir```.

If running with ```--source default```, then provide a ```--sample_sheet``` that contains columns "strain" and "bam_path" (include these headers and absolute paths to bam files). The user can specify an ```--outdir``` where all genomes and statistics will be deposited. 

## --ex_master
If running with ```--ext_master```, and using the master sheet for merging and re-assembling, then provide either a ```--sample_sheet``` or ```--raw_dir``` in addition to the ```--ext_master``` sheet for running in ```--source default``` or ```umd```, respectively. If ```--ext_master``` is run with NO matching strains in raw.dir or the provided sample_sheet, the following error will occur:"ERROR ~ Unknown method invocation `contains` on PoisonPill type -- Did you mean? toString"

This pipeline can also be run in ```--source default``` with a ```--sample_sheet``` that contains bams from multiple runs of the same strain, and merging of sequencing data for the same strain will still occur as the pipeline checks for multiple entries of the same strain.

## --type & --blobtools 
To filter out contigs that might be constructed from non-Nematoda DNA, use parameter ```--blobtools``` set to "yes" when running with ```--type reads```. If the user desires to only run blobtools on already assembled contigous genomes in the FASTA format, the user can specify ```--type assembly``` and exclude the ```--blobtools``` parameter and provide the absolute path to a TSV file with header: strain, asm_fa, species, bam, r_stats. The absolute paths must be provided for asm_fa, bam, and r_stats, and the TSV should be supplied to parameters ```--sample_sheet```. If working with an assembly that was created from merging bams, please provide multiple entries for the same assembly - one entry for each bam file.

When running in Rockfish, use -profile rockfish and it will use the conf/rockfish.config configuration file. Use flag "-resume" to resume an analysis and retrieve any cached data.

## Use cases
### Running in --source umd with --blobtools
```nextflow run main_wBlobToolKit.nf --source default --type reads --sample_sheet --raw_dir 2025_newData --blobtools yes -profile rockfish```

### Running in --source default with no blobtools
```nextflow run main_wBlobToolKit.nf --source default --type reads --sample_sheet /<absolutePath>/strain_bamPaths.tsv --outdir 2025_assemblies -profile rockfish```

### Running in --source default with an --ext_master sheet provided to merge bam files of the same strain, and then running blobtools on the assemblies
```nextflow run main_wBlobToolKit.nf --source default --type reads --sample_sheet /<absolutePath>/strain_bamPaths.tsv --outdir assembleAndBlob --blobtools yes -profile rockfish```

### Running blobtools on already assembled genomes
```nextflow run main_wBlobToolKit.nf --source default --type assembly --sample_sheet /<absolutePath>/strain_asmPaths.tsv --outdir blob_assembly -profile rockfish```
