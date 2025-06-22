
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
 
    // seqrun_file = '"https://docs.google.com/spreadsheets/d/1CpSpzU1p-WtGKIMBK99DL5AeZb-A8QrHPuLkM_fAuEY/export?gid=484600292&format=csv"'
    
    // Order of URLs is CB, CT, CE 
    seqrun_files = Channel.of([
        '"https://docs.google.com/spreadsheets/d/1IJHMLwuaxS_sEO31TyK5NLxPX7_qSd0bHNKverAv8-0/export?gid=1386059628&format=tsv"',
        '"https://docs.google.com/spreadsheets/d/1mqXOlUX7UeiPBe8jfAwFZnqlzhb7X-eKGK_TydT7Gx4/export?gid=1642815395&format=tsv"',
        '"https://docs.google.com/spreadsheets/d/1Rts4CZxkDiid3hux7EpE7CBAQfole6oWQs61dorYBX0/export?gid=538533765&format=tsv"',
        '"https://docs.google.com/spreadsheets/d/1-5EiJHmMqVBm0Emj_wekdDRv_-dAdmmwS9poz_Jxb90/export?gid=578306341&format=tsv"'
    ])
    //seqrun_file_cn = '""'
    //seqrun_file_unnamed = '""'
    seq_ch = get_seqrun(seqrun_files)
            .splitCsv(sep: "\t",header: true)
            .map {row -> [row.sample, row.species] }
            .collect(flat: false)
            .flatten()
            .buffer(size: 2)
            // NEED TO FIX VX34 FOR BRIGGSAE!?!?!

            //.view()
            //.map { rows -> ["PB306", "CE"] + rows } 
            //.concat(Channel.of(["PB306", "CE"]))
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
            // .view { row -> "MASTER: ${row[0]} - ${row[1]}" } // Debugging ext_master	    
        
        // bam_ch.view { row -> "BAM: ${row[0]} - ${row[1]}" } // Debugging bam input
    
        master_join = master_ch
            .groupTuple() 
            .join(bam_ch)
	    .map { row -> [row[0], row[1] + row[2]] }
	    // .view { row -> "JOINED: ${row[0]} - ${row[1]}" } // Debugging join output
 
        master_keys_ch = master_join.map { it[0] }.collect().view() // Collects strain keys

        bam_nmerge = bam_ch.filter { row ->
            def master_keys = master_keys_ch.get() // Ensures `collect()` completes before use
            !master_keys.contains(row[0]) // Filters out matching rows
        }

        // master_keys_ch = master_join.map { it[0] }.collect()
    
        // bam_nmerge = bam_ch.filter { row ->
        // def master_keys = master_keys_ch.get().toSet() // Waits until `collect()` completes
        // !master_keys.contains(row[0]) // Filters out matching rows
        // }

        merged_bams = merge_bam(master_join)
    
        bam_ch_merged = merged_bams.merged
            .mix(bam_nmerge)
	    // .view()
    }
    
   
    mapped_sp_bam = bam_ch_merged
                        .join(seq_ch)
                        //.view()
    
    markdup(mapped_sp_bam)

    rstat_ch = markdup.out.rstat.collectFile(name: "rstat_out.txt")
    
    fastafilt(markdup.out.uniq)

    funiq_ch = fastafilt.out.funiq
                   .filter { it[1].size() > 0 }

    assemble(funiq_ch)
    
    astat_ch = assemble.out.astat.collectFile(name: "astat_out.txt", keepHeader: true, skip: 1)
    
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
    val(seqrun_files)
    
    output:
    path("sp2str_table.tsv"), emit: seqrun

    script:
    """
    wget -O sp.csv ${seqrun_files[0]}
    wget -O sp2.csv ${seqrun_files[1]}
    wget -O sp3.csv ${seqrun_files[2]}
    wget -O sp4.csv ${seqrun_files[3]}

    awk -F'\t' -v OFS='\t' 'NR != 1 {print \$3,"CB"}' sp.csv > cb.tsv
    awk -F'\t' -v OFS='\t' 'NR != 1 {print \$3,"CT"}' sp2.csv > ct.tsv
    awk -F'\t' -v OFS='\t' 'NR != 1 {print \$3,"CE"}' sp3.csv > ce.tsv
    awk -F'\t' -v OFS='\t' 'NR != 1 {print \$3,"CN"}' sp4.csv > cn.tsv

    echo -e "sample\tspecies" > header.tsv
    cat header.tsv cb.tsv ct.tsv ce.tsv cn.tsv | uniq  > sp2str_table.tsv

    #awk -F ',' -v OFS='\t' '{print \$4,\$5}' sp.csv |\
    #sed 's/C\\.[[:space:]]*elegans/CE/' | \
    #sed 's/C\\.[[:space:]]*tropicalis/CT/' | \
    #sed 's/C\\.[[:space:]]*briggsae/CB/' | \
    #sed 's/C\\.[[:space:]]*nigoni/CN/' | \
    #s ed 's/\\.[[:space:]]/\\./g' |
    #awk -F'\t' '{sub(/_.*/,"",\$1); print \$1 "\t" \$2}' | \
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
    echo -e "strain\tbam_path" > header.txt
    
    # Extract strain name and corresponding .bam path from full path in UMD source 
    printf "%s\\n" ${input_dir}/*/*/*.bam | sed 's|${input_dir}/||' | sed 's|/| |' | awk -F' ' -v OFS='\t' '{sub(/_.*/,"",\$1); print \$1}' >> strains.txt
    printf "%s\\n" $input_path/${input_dir}/*/*/*.bam > bams.txt
    paste -d'\t' strains.txt bams.txt > temp1.txt
    cat header.txt temp1.txt > sample_sheet_bam.txt
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
    tuple val(strain), path("$species/assemblies/${uniq.baseName}.${strain}.inbred.asm.bp.p_ctg.fa"), val(species), emit: asm
    path("$species/asm_stat/${uniq.baseName}.${strain}.inbred.asm.bp.p_ctg.fa.stats"), emit: astat

    script:
    """
    mkdir -p ${species}/asm_stat/
    mkdir -p ${species}/assemblies/
    hifiasm -f0 -l0 -t 48 -o ${uniq.baseName}.${strain}.inbred.asm $uniq
    awk '/^S/{print ">"\$2;print \$3}' ${uniq.baseName}.${strain}.inbred.asm.bp.p_ctg.gfa  > $species/assemblies/${uniq.baseName}.${strain}.inbred.asm.bp.p_ctg.fa
    stats.sh -format=6 -in=$species/assemblies/${uniq.baseName}.${strain}.inbred.asm.bp.p_ctg.fa -format=6 -gcformat=0 | awk -v strain=$strain -v OFS='\t' 'NR == 1 {print "strain", \$0} NR > 1 {print strain, \$0}' > $species/asm_stat/${uniq.baseName}.${strain}.inbred.asm.bp.p_ctg.fa.stats
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
    cat $seqflat > ${odir}_SP_CONTENT.txt 
    sort -k1,1 ${odir}_SP_CONTENT.txt | grep -v "(" | awk '\$2 ~ /^(CE|CB|CT|CN)\$/' > ${odir}_sorted.sp2str.txt
    join -t \$'\t' all_body.txt ${odir}_sorted.sp2str.txt > ${odir}_all_body_sp.txt
    cat all_header.txt ${odir}_all_body_sp.txt > ${odir}_all_stats.txt
    """
    
    stub:
    """
    touch ${odir}_all_stats.txt
    """
}
