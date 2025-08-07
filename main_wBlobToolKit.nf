
nextflow.enable.dsl=2

/*
    Parameters
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


/*
    Workflow
*/
workflow {
     
    // Order of URLs is CB, CT, CE 
    seqrun_files = Channel.of([
        '"https://docs.google.com/spreadsheets/d/1IJHMLwuaxS_sEO31TyK5NLxPX7_qSd0bHNKverAv8-0/export?gid=1386059628&format=tsv"',
        '"https://docs.google.com/spreadsheets/d/1mqXOlUX7UeiPBe8jfAwFZnqlzhb7X-eKGK_TydT7Gx4/export?gid=1642815395&format=tsv"',
        '"https://docs.google.com/spreadsheets/d/1Rts4CZxkDiid3hux7EpE7CBAQfole6oWQs61dorYBX0/export?gid=538533765&format=tsv"',
        '"https://docs.google.com/spreadsheets/d/1-5EiJHmMqVBm0Emj_wekdDRv_-dAdmmwS9poz_Jxb90/export?gid=578306341&format=tsv"'
    ])
    // NEED TO ADD A SHEET THAT CONTAINS SP27 AND SP30


    seq_ch = get_seqrun(seqrun_files)                                           
            .splitCsv(sep: "\t",header: true)
            .map {row -> [row.sample, row.species] }
            .collect(flat: false)
            .flatten()
            .buffer(size: 2)

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

        merged_bams = merge_bam(grouped_bam)

        bam_ch_merged = merged_bams
            .mix(single_ch)

    } else { 
        master_ch = Channel.fromPath(params.ext_master, checkIfExists: true)
	    .splitCsv(sep: "\t",header: true)
            .map { row -> [row.strain, row.bam_path] }
            
        master_join = master_ch
            .groupTuple() 
            .join(bam_ch) // this is an inner join - so if there is not a strain to join by, then it is dropped - so this only contains strains that have been previously sequenced
	    .map { row -> [row[0], row[1] + row[2]] } 
 
        master_keys_ch = master_join.map { it[0] }.collect().view() // Collects strain keys - prints out the name of the strains that have been sequenced previously and whose data will be merged with current run

        bam_nmerge = bam_ch.filter { row ->
            def master_keys = master_keys_ch.get() // Ensures `collect()` completes before use
            !master_keys.contains(row[0]) // Filters out matching rows (i.e, strains that have been sequenced multiple times and are in master_join)
        }

        merged_bams = merge_bam(master_join)
    
        bam_ch_merged = merged_bams.merged
            .mix(bam_nmerge) // .mix is like dplyr::bind_rows - binds together the table of strain with bams (some merged, some not merged))
    }
    
   
    mapped_sp_bam = bam_ch_merged
                        .join(seq_ch) // adding species resolution - an inner_join so we must always use species sheets that contain species resolution for all of the strains we are working with
    
    markdup(mapped_sp_bam)                                                       // DO WE REALLY NEED TO PUBLISH THESE FASTAS????????????????????????????????????????????????????????????????????????????????????

    rstat_ch = markdup.out.rstat.collectFile(name: "rstat_out.txt")
    
    fastafilt(markdup.out.uniq)

    funiq_ch = fastafilt.out.funiq
                   .filter { it[1].size() > 0 }

    assemble(funiq_ch)
    
    astat_ch = assemble.out.astat.collectFile(name: "astat_out.txt", keepHeader: true, skip: 1) // keeps the header for the first file, but then appends everythign but the header for the subsequent files (essentially dpyr::bind_rows in R)
    
    seq_flat = seq_ch                                                                                                           /// fix how seq_flat is created!!!!!
                .map { row -> row.join('\t') }
                .map { rows -> ["sp2str.txt"] + rows }
                .collect()
                .collectFile(name: "sp2str.txt", keepHeader: false, newLine: true)
                //.view()
    
    /*
    seq_flat = seq_ch
                .map { row -> ["sp2str.txt", row.join('\t')] }
                .collect()
                .collectFile(name: "sp2str.txt", keepHeader: false, newLine: true)
                */



    if (params.source == "default") {
        gatherstats(rstat_ch, astat_ch, params.outdir, seq_flat) // I NEED TO FIX the line where the second column is grepped for only CE, CT, CN, or CB -  AS SP. 27 AND SP. 30 WILL NOT BE INCLUDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    } else {
        gatherstats(rstat_ch, astat_ch, params.raw_dir, seq_flat)
    }   

    blob_ch = assemble.out.asm
                .join(mapped_sp_bam)          // join by strain
                .map { strain, asm_fa, species1, bam, species2 -> tuple(strain, asm_fa, species1, bam) }  // drop duplicate species2
                .view()
    
    blobtools(blob_ch)                          // REMOVE SECOND MKDIR COMMAND IN PROCESS

    filtasm_ch = blobtools.out.filtasm
                .map { strain, filt_asm, species -> tuple(strain, filt_asm, species) }
                .view()

    busco(filtasm_ch)

    filt_asm_stat_ch = blobtools.out.filtAsmStat.collectFile(name: "filt_astat_out.txt", keepHeader: true, skip: 1)
    busco_out_ch = busco.out.bsco.collectFile(name: "busco_scores.tsv", keepHeader: false, skip: 1)
                    .view()

    gatherstatsFiltered(filt_asm_stat_ch, busco_out_ch, seq_flat, rstat_ch)

    // add a BRAKER3 process - make a different nextflow script that runs BRAKER3, creates proteomes, runs BUSCO on these proteomes, and creates proteomes that are the longest isoform for OrthoFinder
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

    output:
    path("${odir}_all_stats.txt")
    // path("${odir}_SP_CONTENT.txt")
    // path("${odir}_sorted.sp2str.txt")
    // path("all_body.txt")

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


process blobtools {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )

    label 'blobtools'

    input:
    tuple val(strain), path(asm_fa), val(species), path(bam)

    output:
    tuple val(strain), path("${species}/assemblies/filtered/${strain}/${asm_fa.baseName}.filtered.fa"), val(species), emit: filtasm
    path("${species}/asm_stat/filtered/${strain}/${asm_fa.baseName}.filtered.fa.stats"), emit: filtAsmStat


    script:
    """
    mkdir -p ${species}/asm_stat/filtered/${strain}/png
    mkdir -p ${species}/assemblies/filtered/${strain}


    samtools fastq -@ 36 ${bam} | gzip - > ${species}/assemblies/filtered/${strain}/${bam.baseName}.fq.gz
    
    minimap2 -ax map-hifi ${asm_fa} ${species}/assemblies/filtered/${strain}/${bam.baseName}.fq.gz | samtools sort -@ 36 -o ${species}/asm_stat/filtered/${strain}/${bam.baseName}_coverage.bam

    samtools index -c ${species}/asm_stat/filtered/${strain}/${bam.baseName}_coverage.bam


    # Creating a BlobDir:
    blobtools create \
        --fasta ${asm_fa} \
        ${species}/asm_stat/filtered/${strain}_blobDir


    # BLASTing assembly contigs:
    blastn -db /vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/blobtools/core_nt/core_nt \
        -query ${asm_fa} \
        -outfmt "6 qseqid staxids bitscore std" \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads 36 \
        -out ${species}/asm_stat/filtered/${strain}/${strain}_asm_blast.out


    # Adding coverage and BLAST hits to BlobDir:
    blobtools add \
        --hits ${species}/asm_stat/filtered/${strain}/${strain}_asm_blast.out \
        --taxrule bestsumorder \
        --taxdump /vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/blobtools/taxdump \
        --cov ${species}/asm_stat/filtered/${strain}/${bam.baseName}_coverage.bam \
        ${species}/asm_stat/filtered/${strain}_blobDir


    # Filtering out non-Nematoda contigs:                                                             SHOULD I INCLUDE NO-HIT CONTIGS AS WELL???
    blobtools filter \
        --param bestsumorder_phylum--Inv=Nematoda \
        --output ${species}/asm_stat/filtered/${strain}_blobDir/${strain}_nematoda_only_blobDir \
        --fasta ${asm_fa} \
        ${species}/asm_stat/filtered/${strain}_blobDir


    # Creating plots before and after filtering:
    blobtools view \
        --plot \
        --view blob \
        --out ${species}/asm_stat/filtered/${strain}/png \
        ${species}/asm_stat/filtered/${strain}_blobDir

    blobtools view \
        --plot \
        --view blob \
        --out ${species}/asm_stat/filtered/${strain}/png \
        ${species}/asm_stat/filtered/${strain}_blobDir/${strain}_nematoda_only_blobDir
    

    cp ${asm_fa.baseName}.filtered.fa  ${species}/assemblies/filtered/${strain}/${asm_fa.baseName}.filtered.fa

    stats.sh -format=6 -in=${species}/assemblies/filtered/${strain}/${asm_fa.baseName}.filtered.fa -format=6 -gcformat=0 | awk -v strain=$strain -v OFS='\t' 'NR == 1 {print "strain", \$0} NR > 1 {print strain, \$0}' > ${species}/asm_stat/filtered/${strain}/${asm_fa.baseName}.filtered.fa.stats
    """
}


process busco {

   publishDir(
        path: "${params.output}",
        mode: 'copy',
    )
    
    label 'busco'

    input:
    tuple val(strain), path(filt_asm), val(species)

    output:
    path("${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco/${filt_asm.baseName}.busco.stat.tsv"), emit: bsco
    // path("${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco/tmp.tsv")
    // path("${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco/tmp2.tsv")


    script:
    """
    busco -i $filt_asm -c 12 -m genome -l /vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/busco_downloads/lineages/nematoda_odb10/ -o ${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco 

    echo -e "strain\tbusco_completeness\tasm_path" > header.tsv
    grep "C:" ${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco/short_summary.specific.nematoda_odb10.${filt_asm.baseName}.busco.txt > ${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco/tmp.tsv
    awk '{ match(\$0, /C:([0-9.]+)%/, a); print a[1] }' ${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco/tmp.tsv > ${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco/tmp2.tsv 
    paste -d '\t' <(echo "$strain") ${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco/tmp2.tsv <(echo "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/${params.output}/${species}/assemblies/filtered/${strain}/${filt_asm.baseName}.fa") > strain_busco.tsv
    
    cat header.tsv strain_busco.tsv > ${species}/asm_stat/filtered/${strain}/${filt_asm.baseName}.busco/${filt_asm.baseName}.busco.stat.tsv
    """
}


process gatherstatsFiltered {

    publishDir(
        path: "${params.output}",
        mode: 'copy',
    )
    
    label 'gather_statsFiltered'

    input:
    val(filt_asm_stat)
    val(filt_asm_busco)
    val(seqflat)
    val(rstat)

    output:
    path("${params.outdir}_filtered_asm_stats.txt") 

    script:
    """
    cat $filt_asm_stat > filt_asm_stat.txt
    cat $filt_asm_busco > busco_asmPath.txt
    cat $rstat | sort -k1,1 > read_stats.txt
    cat $seqflat | grep -v "sp2str.txt" | awk -F'\t' 'NF && \$1 != "" {print \$0}' > species.txt 
    sort -k1,1 species.txt > species_sorted.txt

    grep "^strain" filt_asm_stat.txt | sed 's/N50/X50/g' | sed 's/L50/N50/g' | sed 's/X50/L50/g' | sed 's/N90/X90/g' | sed 's/L90/N90/g' | sed 's/X90/L90/g' | sed 's|#||' | awk -v OFS='\t' '{print \$0, "yield", "mead_readlen", "species", "genome_busco", "asm_path"}' > header.txt
    grep -v "^strain" filt_asm_stat.txt | sort -k1,1 > body_filt_asm_stats.txt

    join -t \$'\t' body_filt_asm_stats.txt read_stats.txt > body_filt_asm_read_stats.txt
    join -t \$'\t' body_filt_asm_read_stats.txt species_sorted.txt > body_filt_asm_stats_species.txt

    #while IFS=\$'\t' read -r line; do
    #   strain=\$(echo "\$line" | awk '{print \$1}')
    #    busco_asmPath=\$(grep -m1 -w "\$strain" busco_asmPath.txt | awk -v OFS='\t' '{print \$2,\$3}')
    #    echo -e "\${line}\t\${busco_asmPath}" >> final_body.txt
    #done < body_filt_asm_stats_species.txt

    awk 'BEGIN { OFS="\t" }
    NR==FNR { busco[\$1]=\$2"\t"\$3; next }
    { print \$0, busco[\$1] }' busco_asmPath.txt body_filt_asm_stats_species.txt > final_body.txt


    cat header.txt final_body.txt | awk -v OFS='\t' '{print \$23, \$1, \$3, \$5, \$16, \$9, \$10, \$13, \$14, \$17, \$18, \$19, \$21, \$22, \$24, \$25}' > ${params.outdir}_filtered_asm_stats.txt
    """
}

