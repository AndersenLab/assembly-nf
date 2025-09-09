# assembly-nf
Nextflow pipeline for hifiasm *de novo* assembly and sequencing quality control.

By default, the pipeline will assemble genomes in "inbred" mode for hifiasm

## --source
If running with --source umd, then provide the folder in /vast/eande106/data/transfer/raw where sequencing data has been deposited (e.g. --raw_dir 20250314_PacBio). When running in --source umd, the output directory will be set to "${raw_dir}-assembly" even if you specify an --outdir.

If running with --source default, then provide a --sample_sheet that contains columns "strain" and "bam_path" (include these headers and absolute paths to bam files). The user can specify an --outdir where all genomes and statistics will be deposited. 

## --ex_master
If running with --ext_master, and using the master sheet for merging and re-assembling, then provide either a --sample_sheet or --raw_dir in addition to the --ext_master sheet for running in --source default or umd, respectively. If --ext_master is run with NO matching strains in raw.dir or the provided sample_sheet, the following error will occur:"ERROR ~ Unknown method invocation `contains` on PoisonPill type -- Did you mean? toString"

This pipeline can also be run in --source default with a --sample_sheet that contains bams from multiple runs of the same strain, and merging of sequencing data for the same strain will still occur as the pipeline checks for multiple entries of the same strain.

## --type & --blobtools 
To filter out contigs that might be constructed from non-Nematoda DNA, use parameter "--blobtools" set to "yes" when running with '--type reads'. If the user desires to only run blobtools on already assembled contigous genomes in the FASTA format, the user can specify '--type assembly' and exclude the --blobtools parameter and provide the absolute path to a TSV file with header: strain, asm_fa, species, bam, r_stats. The absolute paths must be provided for asm_fa, bam, and r_stats, and the TSV should be supplied to parameters --sample_sheet. If working with an assembly that was created from merging bams, please provide multiple entries for the same assembly - one entry for each bam file.

When running in Rockfish, use -profile rockfish and it will use the conf/rockfish.config configuration file. Use flag "-resume" to resume an analysis and retrieve any cached data.

## Examples
### Running in --source umd 

### Running in --source default

### Running in --source default with an --ext_master sheet proided to merge bam files of the same strain, and then running blobtools on the assemblies

### Running blobtools on already assembled genomes
