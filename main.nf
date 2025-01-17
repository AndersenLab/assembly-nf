nextflow.enable.dsl=2

/*
    Params
*/

date = new Date().format('yyyyMMdd')

if (params.debug) {
    println """
    *** Using debug mode ***
    """
    params.output = "assembly-${date}-debug"
    params.pbdata = "${workflow.projectDir}/test_data"
} else {
    params.output = "assembly-${date}"
    params.pbdata = "${params.data_path}/${params.raw_dir}"
}

//log.info("${params.data_path}")

def log_summary() {
    // Corrected log summary function to print information instead of recursive call
    log.info("Workflow summary: \n" + 
             "Debug mode: ${params.debug}\n" + 
             "Output directory: ${params.output}\n" + 
             "PB Data path: ${params.pbdata}")
    
    // Show help if requested
    if (params.help) {
        log.info("Help requested, exiting.")
        exit 1
    }
}

// Start the workflow
workflow {
    // Call the gensheet process and store the result in a variable
    input_dir = file(params.pbdata)

    file_list = gensheet(input_dir,params.data_path)
    
    bam_ch = gensheet.out.bam
        .splitCsv(sep: "\t",header: true)
	.view { row -> "${row.strain} - ${row.bam_path}" }

    markdup(bam_ch)
    // hifiasm(markdup.out.uniq))
}

process gensheet {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )

    label 'create_sample_sheet'

    input:
    path input_dir
    val input_path

    output:
    path("sample_sheet_bam.txt"), emit: bam

    script:
    """
    echo -e "strain\tbam_path" > sample_sheet_bam.txt
    # Ensure the input directory is absolute
    # input_dir=\$(realpath ${input_dir})
    # echo "$input_path"
    # List directories and strip paths
    printf "%s\\n" ${input_dir}/*/ | sed 's/ /\\\\n/g' | sed 's|.*/\\([^/]*\\)/[^/]*|\\1|' > strains.tmp

    #printf "%s\\n" ${input_dir}/*/*/*/*.bam | sed 's/ /\\\\n/g' > bams.tmp
    # List .bam files and replace spaces with newlines
    printf "%s\\n" ${input_dir}/*/*/*.bam | sed 's/ /\\\\n/g' | awk -v append_path="${input_path}" -v OFS="/" '{print append_path,\$0}' > bams.tmp

    paste strains.tmp bams.tmp >> sample_sheet_bam.txt
    """
}

process markdup {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )

    label 'pb_mark_duplicates'

    input:
    tuple val(strain), path(bam) 

    output:
    tuple val(strain), path("${bam.baseName}.uniq.fasta"), emit: uniq

    script:
    """
    pbmarkdup $bam ${bam.baseName}.uniq.fasta --dup-file ${bam.baseName}.dups.fasta
    """
}
