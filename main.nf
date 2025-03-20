
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
    
    seqrun_file = '"https://docs.google.com/spreadsheets/d/1CpSpzU1p-WtGKIMBK99DL5AeZb-A8QrHPuLkM_fAuEY/export?gid=484600292&format=csv"'
    seq_ch = get_seqrun(seqrun_file)
            .splitCsv(sep: "\t",header: true)
            .map {row -> [row.sample, row.species] }
            .collect(flat: false)
            .flatten()
            .buffer(size: 2)
            //.map { rows -> ["PB306", "CE"] + rows }
            .concat(Channel.of(["PB306", "CE"]))
            //.view()
            //.map { row -> [row[0], row[1]] }
            //.view{ row -> "$row[0] - $row[2]}

    // Call the gensheet process and store the result in a variable
    if (params.source == "umd") {

        input_dir = file(params.pbdata)
        file_list = gensheet(input_dir,params.data_path)
    
        bam_ch = gensheet.out.bam
            .splitCsv(sep: "\t",header: true)
	    .map { row -> [row.strain, row.bam_path] }

    } else {
          
         bam_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)
            .ifEmpty { exit 1, "sample sheet not found" }
            .splitCsv(sep: "\t",header: true)
            .map { row -> [row.strain, row.bam_path] }
    }
    
    if (params.ext_master == null) {
        grouped_bam = bam_ch
            .groupTuple()
            .filter { row -> row[1].size() > 1 }
        
        single_ch = bam_ch
                .groupTuple()
                .filter { row -> row[1].size() == 1 }
                .map { row -> row[1] = row[1].first()
                              return row }  
                //.view()

        merged_bams = merge_bam(grouped_bam)

        bam_ch_merged = merged_bams
            .mix(single_ch)
             //.view()

    } else { 
        master_ch = Channel.fromPath(params.ext_master, checkIfExists: true)
	    .splitCsv(sep: "\t",header: true)
            .map { row -> [row.strain, row.bam_path] }
	    //.view { row -> "${row.strain} - ${row.bam_path} - ${row.run}" }

        master_join = master_ch
            .groupTuple() 
            .join(bam_ch)
	    .map { row -> [row[0], row[1] + row[2]] }
	    //.view()
 
        master_keys_ch = master_join.map { it[0] }.collect()

    
        bam_nmerge = bam_ch.filter { row ->
        def master_keys = master_keys_ch.get().toSet() // Waits until `collect()` completes
        !master_keys.contains(row[0]) // Filters out matching rows
        }

        merged_bams = merge_bam(master_join)
    
        bam_ch_merged = merged_bams.merged
            .mix(bam_nmerge)
	    //.view()
    }
    
    //seq_ch.view()
    mapped_sp_bam = bam_ch_merged
                        //.view()
                        .join(seq_ch)
                        //.view()
    
    markdup(mapped_sp_bam)

    rstat_ch = markdup.out.rstat.collectFile(name: "rstat_out.txt")
    
    fastafilt(markdup.out.uniq)

    funiq_ch = fastafilt.out.funiq
                   .filter { it[1].size() > 0 }

    assemble(funiq_ch)
    
    astat_ch = assemble.out.astat.collectFile(name: "astat_out.txt", keepHeader: true, skip: 1)
    
    // seq_ch.view()
    seq_flat = seq_ch
                .map { row -> row.join('\t') }
                //.map { row -> "${row[0]}\t${row[1]}" }
                //.view()
                //.map { it.join(",") }
                //.view()
                //.collect()
                .map { rows -> ["sp2str.txt"] + rows }
                .collect()
                .collectFile(name: "sp2str.txt", keepHeader: false, newLine: true)
                //.view()

    if (params.source == "default") {
        gatherstats(rstat_ch, astat_ch, params.outdir, seq_flat)
    } else {
        gatherstats(rstat_ch, astat_ch, params.raw_dir, seq_flat)
    }
}

process get_seqrun {
    
    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )

    input:
    val(seqrun)
    
    output:
    path("sp2str_table.tsv"), emit: seqrun

    script:
    """
    #awk '{print \$0}' $seqrun > sp2str_table.tsv
    wget -O sp.csv $seqrun
    awk -F ',' -v OFS='\t' '{print \$4,\$5}' sp.csv |\
    sed 's/C\\.[[:space:]]*elegans/CE/' | \
    sed 's/C\\.[[:space:]]*tropicalis/CT/' | \
    sed 's/C\\.[[:space:]]*briggsae/CB/' | \
    sed 's/C\\.[[:space:]]*nigoni/CN/' | \
    sed 's/\\.[[:space:]]/\\./g' |
    awk -F'\t' '{sub(/_.*/,"",\$1); print \$1 "\t" \$2}' | \
    uniq  > sp2str_table.tsv
    """
}

process gensheet {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )

    label 'create_sample_sheet'

    input:
    path(input_dir)
    val(input_path)

    output:
    path("sample_sheet_bam.txt"), emit: bam

    script:
    """
    echo -e "strain\tbam_path" > sample_sheet_bam.txt
    # Ensure the input directory is absolute
    # input_dir=\$(realpath ${input_dir})
    # echo "$input_path"
    # List directories and strip paths
    printf "%s\\n" ${input_dir}/*/ | sed 's/ /\\\\n/g' | sed 's|.*/\\([^/]*\\)/[^/]*|\\1|' | sed 's/_.*//' > strains.tmp  ############ fix to extract strain from bam path 

    #printf "%s\\n" ${input_dir}/*/*/*/*.bam | sed 's/ /\\\\n/g' > bams.tmp
    # List .bam files and replace spaces with newlines
    printf "%s\\n" ${input_dir}/*/*/*.bam | sed 's/ /\\\\n/g' | awk -v append_path="${input_path}" -v OFS="/" '{print append_path,\$0}' > bams.tmp

    paste strains.tmp bams.tmp >> sample_sheet_bam.txt
    """
}

process merge_bam {

    label 'merge'
    executor workflow.stubRun ? 'local' : 'slurm'

    input:
    tuple val(strain), path(bam)

    output:
    tuple val(strain), file("new_${strain}.bam"), emit: merged

    script:
    """
    samtools merge -o new_${strain}.bam $bam
    #sambamba index --nthreads=${task.cpus} new_${strain}.bam
    """

    stub:
    """
    touch new_${strain}.bam
    """
}

process markdup {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )

    label 'pb_mark_duplicates'

    input:
    tuple val(strain), path(bam), val(species)

    output:
    tuple val(strain), path("${bam.baseName}.uniq.fasta"), val(species), emit: uniq
    path("${species}/read_stat/${bam.baseName}.uniq.fasta.read_stats.txt"), emit: rstat

    script:
    """
    mkdir -p ${species}/read_stat/
    pbmarkdup $bam ${bam.baseName}.uniq.fasta --dup-file ${bam.baseName}.dups.fasta
    count="\$(grep '^>' ${bam.baseName}.uniq.fasta | wc -l)"; echo \$count > read_yield.txt
    grep -v "^>" ${bam.baseName}.uniq.fasta | awk '{total+=length(\$0); count++} END {if(count>0) print total/count; else print "No sequences found"}' > read_avglen.txt
    paste -d '\t' read_yield.txt read_avglen.txt | awk -v strain=$strain -v OFS='\t' '{print strain,\$0}' > ${species}/read_stat/${bam.baseName}.uniq.fasta.read_stats.txt
    #grep "^>" ${bam.baseName}.uniq.fasta > ${bam.baseName}.read_count.txt
    """
    
    stub:
    """
    touch ${bam.baseName}.uniq.fasta
    touch ${bam.baseName}.uniq.fasta.read_stats.txt
    """

}


process fastafilt {

    input:
    tuple val(strain), path(uniq), val(species)  // Input: strain and path to FASTA file

    output:
    tuple val(strain), path("${uniq.baseName}.filtered.fasta"), val(species), emit: funiq
    
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
    
    stub:
    """
    touch ${uniq.baseName}.filtered.fasta
    """

}

process assemble {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )

    label 'pb_assemble'

    input:
    tuple val(strain), path(uniq), val(species)
    
    output:
    tuple val(strain), path("$species/assemblies/${uniq.baseName}.inbred.asm.bp.p_ctg.fa"), val(species), emit: asm
    path("$species/asm_stat/${uniq.baseName}.inbred.asm.bp.p_ctg.fa.stats"), emit: astat

    script:
    """
    mkdir -p ${species}/asm_stat/
    mkdir -p ${species}/assemblies/
    hifiasm -f0 -l0 -t 48 -o ${uniq.baseName}.inbred.asm $uniq
    awk '/^S/{print ">"\$2;print \$3}' ${uniq.baseName}.inbred.asm.bp.p_ctg.gfa  > $species/assemblies/${uniq.baseName}.inbred.asm.bp.p_ctg.fa
    stats.sh -format=6 -in=$species/assemblies/${uniq.baseName}.inbred.asm.bp.p_ctg.fa -format=6 -gcformat=0 | awk -v strain=$strain -v OFS='\t' 'NR == 1 {print "strain", \$0} NR > 1 {print strain, \$0}' > $species/asm_stat/${uniq.baseName}.inbred.asm.bp.p_ctg.fa.stats ############################ ADD STRAIN NAME TO THE FINAL ASSEMBLY
    """

    stub:
    """
    touch ${uniq.baseName}.inbred.asm.bp.p_ctg.fa
    touch ${uniq.baseName}.inbred.asm.bp.p_ctg.fa.stats
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
    path(seqflat)
    //tuple val(sp), val(strain)

    output:
    path("${odir}_all_stats.txt")
    path("${odir}_SP_CONTENT.txt")
    path("${odir}_sorted.sp2str.txt")
    path("all_body.txt")

    script:
    """
    cat $astat > assembly_stats.txt
    cat $rstat | sort -k1,1 > read_stats.txt
    grep "^strain" assembly_stats.txt | sed 's/N50/X50/g' | sed 's/L50/N50/g' | sed 's/X50/L50/g' > header_asm.txt
    echo -e "yield\tmean_readlen\tspecies" > header_read.txt
    paste -d '\t' header_asm.txt header_read.txt > all_header.txt
    grep -v "^strain" assembly_stats.txt | sort -k1,1 > body_asm.txt
    join -t \$'\t' body_asm.txt read_stats.txt > all_body.txt
    cat $seqflat > ${odir}_SP_CONTENT.txt ####################### FIX SP 27 AND SP 30 ##############
    sort -k1,1 ${odir}_SP_CONTENT.txt | grep -v "(" | awk '\$2 ~ /^(CE|CB|CT|CN)\$/' > ${odir}_sorted.sp2str.txt
    join -t \$'\t' all_body.txt ${odir}_sorted.sp2str.txt > ${odir}_all_body_sp.txt
    cat all_header.txt ${odir}_all_body_sp.txt > ${odir}_all_stats.txt
    """
    
    stub:
    """
    touch ${odir}_all_stats.txt
    """
}
