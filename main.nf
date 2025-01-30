nextflow.enable.dsl=2

/*
    Params
*/

date = new Date().format('yyyyMMdd')
log.info("Source: ${params.source}")

if (params.debug) {
    println """
    *** Using debug mode ***
    """
    params.output = "assembly-${date}-debug"
    params.pbdata = "${workflow.projectDir}/test_data"
} else {
    if (params.source == 'default') {
        println """
            Running on default source.
            """
        if (params.sample_sheet == null) {
            println """
            Please specify a sample sheet with option --sample_sheet.
            """
            exit 1
        }
        if (params.outdir == null) {
            println """
            Please specify an output directory name with option --outdir.
            """
            exit 1
        }
        params.output = "${params.outdir}-assembly"
        // params.pbdata = "${params.data_path}"
    } else if (params.source == 'umd') {
        println """
            Running on UMD source.
            """
        params.output = "${params.raw_dir}-assembly"
        params.pbdata = "${params.data_path}/${params.raw_dir}"
    } else {
	log.info("Missing source, exiting.")
        exit 1
    }
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
    
    if (params.source == "umd") {
        input_dir = file(params.pbdata)
        file_list = gensheet(input_dir,params.data_path)
    
        bam_ch = gensheet.out.bam
            .splitCsv(sep: "\t",header: true)
	    .view { row -> "${row.strain} - ${row.bam_path}" }
    } else {
         bam_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)
            .ifEmpty { exit 1, "sample sheet not found" }
            .splitCsv(sep: "\t",header: true)
            .view { row -> "${row.strain} - ${row.bam_path}" }
    }

    markdup(bam_ch)

    rstat_ch = markdup.out.rstat.collectFile(name: "rstat_out.txt").view()

    fastafilt(markdup.out.uniq)

    funiq_ch = fastafilt.out.funiq
                   .filter { it[1].size() > 0 }

    assemble(funiq_ch)
    
    // asm_ch, astat_ch = assemble(uniq_ch)
    // uniq_ch, rstat_ch = markdup(bam_ch)
    // rstat_ch = markdup.out.rstat
    astat_ch = assemble.out.astat.collectFile(name: "astat_out.txt", keepHeader: true, skip: 1).view()
    
    if (params.source == "default") {
        gatherstats(rstat_ch, astat_ch, params.outdir)
    } else {
        gatherstats(rstat_ch, astat_ch, params.raw_dir)
    }

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
    #grep "^>" ${bam.baseName}.uniq.fasta > ${bam.baseName}.read_count.txt
    """
    //read_count = file("${params.output}/read_count.txt").text.trim()
}


process fastafilt {

    input:
    tuple val(strain), path(uniq)  // Input: strain and path to FASTA file

    output:
    tuple val(strain), path("${uniq.baseName}.filtered.fasta"), emit: funiq
    
    script:
    """
    # Count the number of reads in the FASTA file using grep
    read_count=\$(grep -c '^>' ${uniq} || echo 0)

    # Check if the read count is greater than 500
    if [ \$read_count -gt 500 ]; then
        cp ${uniq} ${uniq.baseName}.filtered.fasta
    else
        echo "FASTA file '${uniq}' has fewer than 500 reads. Skipping."
        touch ${uniq.baseName}.filtered.fasta
    fi
    """
}

//process fastafilt {
//
//    input:
//    tuple val(strain), path(uniq)  // Input tuple: strain and path to FASTA file
//
//    output:
//    tuple val(strain), path("${uniq.baseName}.filtered.fasta"), emit: funiq  // Output tuple: strain and path to FASTA file if read count > 500
//    tuple val(strain), path('dummy.fasta'), emit: dummy  // Dummy output for cases that don't meet the read count condition
//
//    script:
//    """
//    # Count the number of reads in the FASTA file using grep
//    read_count=\$(grep -c '^>' ${uniq})
//
//    # Check if the read count is greater than 500
//    if [ \$read_count -gt 500 ]; then
//        # If the read count is greater than 500, emit the FASTA file
//        echo "FASTA file '${uniq}' has more than 500 reads. Emitting."
//        cp ${uniq} ${uniq.baseName}.filtered.fasta
//        # Emit the FASTA file that passed the read count check
//        echo "${strain}, ${uniq.baseName}.filtered.fasta"
//    else
//        # If the read count is less than or equal to 500, emit a dummy file or null
//        echo "FASTA file '${uniq}' has fewer than 500 reads. Skipping."
//        touch dummy.fasta
//        # Emit a dummy file to satisfy Nextflow's output requirement
//        echo "${strain}, dummy.fasta"
//    fi
//    """
//}


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
    val(odir)

    output:
    path("${odir}_all_stats.txt")

    script:
    """
    cat $astat > assembly_stats.txt
    cat $rstat | sort -k1,1 > read_stats.txt
    grep "^strain" assembly_stats.txt | sed 's/N50/X50/g' | sed 's/L50/N50/g' | sed 's/X50/L50/g' > header_asm.txt
    echo -e "yield\tmean_readlen" > header_read.txt
    paste -d '\t' header_asm.txt header_read.txt > all_header.txt
    grep -v "^strain" assembly_stats.txt | sort -k1,1 > body_asm.txt
    join -t \$'\t' body_asm.txt read_stats.txt > all_body.txt
    cat all_header.txt all_body.txt > ${odir}_all_stats.txt
    """
}
