include { get_software; q2p; print_error } from './models/utils'
software = get_software(params.software_dir, params.conda_path)
// process guppy {
//     storeDir "${params.outDir}/${params.analysisName}/${sampleName}/01.processFastq"
//     label "assembly_guppy"

// }
nextflow.enable.dsl = 2

process kraken2 {
    storeDir        "${params.out_dir}/${params.analysis_name}/${sample_name}/01.process_fastq"
    label           "assembly_kraken2"
    tag             "${sample_name}"
    debug           true

    input:
        tuple   val(sample_name),
                path(input_fastq)
    output:
        tuple   val(sample_name),
                path("${sample_name}_reads_classified.tsv"),
                path("${sample_name}_classified.txt")

    script:
    """
    export KRAKEN2_DB_PATH=${params.kraken_database};
    ${software['kraken2']} --threads ${params.threads} --output ${sample_name}_reads_classified.tsv --report ${sample_name}_classified.txt --db standard ${input_fastq};
    """

    stub:
    """ touch ${sample_name}_reads_classified.tsv && \
    touch ${sample_name}_classified.txt """
}

process get_most_g {
    beforeScript    "bash /home/a/t.sh"
    storeDir        "${params.out_dir}/${params.analysis_name}/${sample_name}/01.process_fastq"
    label           "assembly_get_most_g"
    tag             "${sample_name}"
    debug           true

    input:
        tuple   val(sample_name),
                path(reads_classified_fp),
                path(classified_fp),
                path(input_fastq),
                val(min_read_len),
                val(min_read_qual)
    output:
        tuple   val(sample_name),
                path("${sample_name}_clean.fastq.gz"),
                path("${sample_name}_tmp.fastq"),
                path("${sample_name}_reads.txt")
    script:
    """
    ${software['python']} ${projectDir}/scripts/get_most_g.py ${reads_classified_fp} ${classified_fp} ${sample_name}_reads.txt;
    awk '{print \$1}' ${sample_name}_reads.txt | ${software['seqtk']} subseq ${input_fastq} -  > ${sample_name}_tmp.fastq;
    ${software['filtlong']} --min_length ${min_read_len} --min_mean_q ${min_read_qual} ${sample_name}_tmp.fastq | gzip -c > ${sample_name}_clean.fastq.gz
    """

    stub:
    """ touch ${sample_name}_clean.fastq.gz &&  \
    touch ${sample_name}_tmp.fastq && \
    touch ${sample_name}_reads.txt """
}

process denovo {
    storeDir    "${params.out_dir}/${params.analysis_name}/${sample_name}/02.assembly"
    label       "assembly_assembly"
    conda       "${params.conda_path}/flye"
    tag         "${sample_name}"
    debug       true

    input:
        tuple   val(sample_name),
                path(clean_fastq)
    output:
        tuple   val(sample_name),
                path("${sample_name}_assembly"),
                path("${sample_name}_assembly/${sample_name}_draft.fasta"),
                path("${sample_name}_assembly/${sample_name}_flye_stat.tsv")
    script:
    """
    flye --nano-raw ${clean_fastq} --out-dir ${sample_name}_assembly --threads ${params.threads}; 
    # mv ./outdir/assembly.fasta ${sample_name}_draft.fasta;
    # mv ./outdir/assembly.fasta ${sample_name}_flye_stat.tsv;
    mv ./${sample_name}_assembly/assembly.fasta ./${sample_name}_assembly/${sample_name}_draft.fasta;
    mv ./${sample_name}_assembly/assembly_info.txt ./${sample_name}_assembly/${sample_name}_flye_stat.tsv;
    """

    stub:
    """ mkdir ${sample_name}_assembly &&  \
    touch ${sample_name}_assembly/${sample_name}_draft.fasta && \
    touch ${sample_name}_assembly/${sample_name}_flye_stat.tsv """
}

process denovo2 {
    // for hybrid assembly
    storeDir    "${params.out_dir}/${params.analysis_name}/${sample_name}/02.assembly"
    label       "assembly_hybrida_assembly"
    conda       "${params.condas_path}/unicycler"
    tag         "${sample_name}"
    
    input:
        tuple   val(sample_name),
                path("ngs_reads1.fastq.gz"),
                // path("ngs_reads2.fastq.gz")
                path("ont.fastq.gz")
    
    output:
        tuple   val(sample_name),
                val("./hybrid_assembly/${sample_name}.fasta")
    
    script:
    assert params.reads2
    def reads2 = params.reads2 ? "-2 ${params.reads2}" : ""
    """
    unicycler -1 ngs_reads1.fastq.gz  ${reads2} -l ont.fastq.gz -o hybrid_assembly -t ${params.threads}
    """

}

process polish{
    storeDir        "${params.out_dir}/${params.analysis_name}/${sample_name}/02.assembly"
    label           "assembly_polish"
    conda           "${params.conda_path}/medaka"
    tag             "${sample_name}"

    input:
        tuple       val(iter),  // how many times medaka will be used to polish the draft reference
                    val(sample_name), 
                    path(clean_fastq), 
                    path(draft_assembly)    // the draft assembly from flye or the input reference by users
        
    output:
        tuple       val(sample_name),
                    path("${sample_name}.reads_map_assembly.${iter}.bam"), // the last loop bam will be reserved and uesd to plot the depth
                    path("${sample_name}.reads_map_assembly.${iter}.bam.bai"),
                    path("${sample_name}.reads_map_assembly.${iter}.hdf"),
                    path("${sample_name}.fasta") // the final fasta used for the subsequent analysis
    
    script:
    """
    draft_fasta=${draft_assembly}
    for i in `seq 1 ${iter}`; do 
        mini_align -i ${clean_fastq} -r \${draft_fasta} -p ${sample_name}.reads_map_assembly.\${i} -t ${params.threads};
        medaka consensus ${sample_name}.reads_map_assembly.\${i}.bam ${sample_name}.reads_map_assembly.\${i}.hdf --threads 2;
        medaka stitch --threads 2 ${sample_name}.reads_map_assembly.\${i}.hdf  \${draft_fasta}  ${sample_name}.\${i}.fasta;
        draft_fasta=${sample_name}.\${i}.fasta
    done

    mv ${sample_name}.${iter}.fasta ${sample_name}.fasta
    """

    stub:
    """touch ${sample_name}.reads_map_assembly.${iter}.hdf && \
    touch ${sample_name}.fasta && \
    touch ${sample_name}.reads_map_assembly.${iter}.bam && \
    touch ${sample_name}.reads_map_assembly.${iter}.bam.bai """
}

/*
process mosdepth {
    storeDir        "${params.out_dir}/${params.analysis_name}/${sample_name}/04.depth"
    label           "assembly_depth"
    conda           "${params.conda_path}/ont"
    tag             "${sample_name}"

    input:
        tuple       val(sample_name), path('sample.bam'), path('bam.bai')
    output:
        tuple       val(sample_name), path()
}
*/

process prokka {
    storeDir        "${params.out_dir}/${params.analysis_name}/${sample_name}/03.annotate"
    label           "assembly_prokka"
    conda           "${params.conda_path}/prokka"
    tag             "${sample_name}"

    input:
        tuple       val(sample_name),
                    path('assembly.fasta')
    
    output:
        tuple       val(sample_name),
                    path("${sample_name}.gff"),
                    path("${sample_name}.pro.fasta"),
                    path("${sample_name}.nuc.fasta")
    
    script:
    def prokka_opts = params.prokka_opts ? params.prokka_opts : ""
    """
    prokka ${prokka_opts} --cpus ${params.threads} --prefix ${sample_name} --outdir ./${sample_name}_prokka assembly.fasta
    mv ./${sample_name}_prokka/${sample_name}.gff ${sample_name}.gff
    mv ./${sample_name}_prokka/${sample_name}.faa ${sample_name}.pro.fasta
    mv ./${sample_name}_prokka/${sample_name}.ffn ${sample_name}.nuc.fasta
    """

    stub:
    """touch ${sample_name}.gff ${sample_name}.pro.fasta ${sample_name}.nuc.fasta"""
}

process blast_pubmlst {
    storeDir        "${params.out_dir}/${params.analysis_name}/${sample_name}/03.annotate"
    label           "assembly_blast_mlst"
    tag             "${sample_name}"

    input:
        tuple       val(sample_name),
                    path('assembly.fasta')
        
    output:
        tuple       val(sample_name),
                    path("${sample_name}_blast_pubmlst.tsv")

    script:
    """
    echo ${projectDir}
    ${software['blastn']}   -query assembly.fasta \
                            -db ${projectDir}/database/pubmlst/pubmlst_alleles.fasta \
                            -max_target_seqs 10000 \
                            -outfmt \"6 qseqid  sseqid slen qstart qend  sstart send  length  pident sstrand\" \
                            > ${sample_name}_blast_pubmlst.tsv
    """

    stub:
    """touch ${sample_name}_blast_pubmlst.tsv"""
}

process get_mlst {
    storeDir        "${params.out_dir}/${params.analysis_name}/${sample_name}/03.annotate"
    label           "assembly_get_mlst"
    tag             "${sample_name}"

    input:
        tuple       val(sample_name),
                    path("blast_out.tsv")
    
    output:
        tuple       val(sample_name),
                    path("${sample_name}_mlst_result.tsv")
    
    script:
    """
    ${software['python']}  ${projectDir}/scripts/get_st_from_blast.py  blast_out.tsv  ${projectDir}/database/pubmlst/st_alleles.json ${sample_name}_mlst_result.tsv
    """

    stub:
    """touch ${sample_name}_mlst_result.tsv"""
}

process annotate {
    storeDir        "${params.out_dir}/${params.analysis_name}/${sample_name}/03.annotate"
    label           "assembly_annotate"
    conda           "${params.conda_path}/abricate"
    tag             "${sample_name}"   

    input:
        tuple       val(sample_name),
                    path('assembly.fasta')
    
    output:
        tuple       val(sample_name),
                    path("${sample_name}_card_result.tsv"),
                    path("${sample_name}_vfdb_result.tsv")
    
    script:
    """
    abricate --datadir ${projectDir}/database --db card assembly.fasta > ${sample_name}_card_result.tsv
    abricate --datadir ${projectDir}/database --db vfdb assembly.fasta > ${sample_name}_vfdb_result.tsv
    """

    stub:
    """touch ${sample_name}_card_result.tsv ${sample_name}_vfdb_result.tsv"""
    
}


workflow {
    // def a = ["barcode17", params.input_fastq]
    // kraken2_out = kraken2(a)
    // kraken2_out.view()
    // kraken2_out.map {
    //     it -> it + [params.input_fastq, params.min_read_len, q2p(params.min_read_qual)]
    // }.set { get_most_g_in }
    // get_most_g_out = get_most_g(get_most_g_in)

    // get_most_g_out.map {
    //     it -> it[0, 1]
    // }.set { denovo_in }

    def denovo_in = ['barcode03', params.input_fastq]

    denovo_out = denovo(denovo_in)  // denovo_out: [name, outdir, draft.fasta, stat.tsv]
    denovo_out.map {
        it -> [params.iter, it[0], params.input_fastq, it[2]]
    }.set { polish_in }

    // polish_in: [polish_number, name, input_fastq, draft.fasta]
    polish_out = polish(polish_in)
    
    polish_out.map {
        it -> [it[0], it[-1]]
    }.set { prokka_in }

    // prokka_in: [name, assembly.fasta]
    // blast_pubmlst_in: [name, assembly.fasta]
    // annotate_in: [name, assembly.fasta]
    annotate_in = blast_pubmlst_in = prokka_in
    // prokka_out: [name, gff, pro.fa, nuc.fa]
    prokka_out = prokka(prokka_in)
    annotate_out = annotate(prokka_in)
    blast_pubmlst_out = blast_pubmlst(blast_pubmlst_in)
    // blast_pubmlst_out: [name, pubmlst_blast.tsv]

    get_mlst_out = get_mlst(blast_pubmlst_out)
    get_mlst_out.view()
    annotate_out.view()

}
