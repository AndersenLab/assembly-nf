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
    rstat_ch = markdup.out.rstat.collectFile(name: "rstat_out.txt").view()

    assemble(markdup.out.uniq)
    
    // asm_ch, astat_ch = assemble(uniq_ch)
    // uniq_ch, rstat_ch = markdup(bam_ch)
    // rstat_ch = markdup.out.rstat
    astat_ch = assemble.out.astat.collectFile(name: "astat_out.txt", keepHeader: true, skip: 1).view()
    gatherstats(rstat_ch, astat_ch)

    // gatherstats(markdup.out.rstat, assemble.out.astat)
    // getstats(assemble.out.asm)

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
    path("${bam.baseName}.uniq.fasta.read_stats.txt"), emit: rstat

    script:
    """
    pbmarkdup $bam ${bam.baseName}.uniq.fasta --dup-file ${bam.baseName}.dups.fasta
    count="\$(grep '^>' ${bam.baseName}.uniq.fasta | wc -l)"; echo \$count > read_yield.txt
    grep -v "^>" ${bam.baseName}.uniq.fasta | awk '{total+=length(\$0); count++} END {if(count>0) print total/count; else print "No sequences found"}' > read_avglen.txt
    paste -d '\t' read_yield.txt read_avglen.txt | awk -v strain=$strain -v OFS='\t' '{print strain,\$0}' > ${bam.baseName}.uniq.fasta.read_stats.txt
    """
}

process assemble {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )

    label 'pb_assemble'

    input:
    tuple val(strain), path(uniq)
    
    output:
    tuple val(strain), path("${uniq.baseName}.inbred.asm.bp.p_ctg.fa"), emit: asm
    path("${uniq.baseName}.inbred.asm.bp.p_ctg.fa.stats"), emit: astat

    script:
    """
    hifiasm -f0 -l0 -t 12 -o ${uniq.baseName}.inbred.asm $uniq
    awk '/^S/{print ">"\$2;print \$3}' ${uniq.baseName}.inbred.asm.bp.p_ctg.gfa > ${uniq.baseName}.inbred.asm.bp.p_ctg.fa
    stats.sh -format=6 -in=${uniq.baseName}.inbred.asm.bp.p_ctg.fa -format=6 -gcformat=0 | awk -v strain=$strain -v OFS='\t' 'NR == 1 {print "strain", \$0} NR > 1 {print strain, \$0}' > ${uniq.baseName}.inbred.asm.bp.p_ctg.fa.stats
    """
}

process gatherstats {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )

    label 'gather_stats'

    input:
    val(rstat)
    val(astat)

    output:
    path("all_stats.txt")

    script:
    """
    cat $astat > assembly_stats.txt
    cat $rstat | sort -k1,1 > read_stats.txt
    grep "^strain" assembly_stats.txt > header_asm.txt
    echo -e "yield\tmean_readlen" > header_read.txt
    paste -d '\t' header_asm.txt header_read.txt > all_header.txt
    grep -v "^strain" assembly_stats.txt | sort -k1,1 > body_asm.txt
    join -t \$'\t' body_asm.txt read_stats.txt > all_body.txt
    cat all_header.txt all_body.txt > all_stats.txt
    """
}
