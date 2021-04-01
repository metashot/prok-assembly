nextflow.enable.dsl=2

process spades {
    tag "${id}"

    publishDir "${params.outdir}/spades" , mode: 'copy' ,
        pattern: "${id}/*"

    publishDir "${params.outdir}/scaffolds" , mode: 'copy' ,
        pattern: "${id}.fa"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.fa"), emit: scaffolds
    path "${id}/*"

    script:
    task_memory_GB = task.memory.toGiga()
    param_spades_k = params.spades_k == 'default' ? "" :  "-k ${params.spades_k}"
    input = params.single_end ? "-s \"$reads\"" : "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    spades.py \
        $input \
        --isolate \
        ${param_spades_k} \
        --threads ${task.cpus} \
        --memory ${task_memory_GB} \
        -o ${id}
    cp ${id}/scaffolds.fasta ${id}.fa
    """
}

process plasmidspades {
    tag "${id}"

    publishDir "${params.outdir}/plasmidspades" , mode: 'copy' ,
        pattern: "${id}/*"

    publishDir "${params.outdir}/scaffolds_plasmids" , mode: 'copy' ,
        pattern: "${id}.fa"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.fa"), emit: scaffolds
    path "${id}/*"

    script:
    task_memory_GB = task.memory.toGiga()
    param_spades_k = params.spades_k == 'default' ? "" :  "-k ${params.spades_k}"
    input = params.single_end ? "-s \"$reads\"" : "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    spades.py \
        $input \
        --plasmid \
        ${param_spades_k} \
        --threads ${task.cpus} \
        --memory ${task_memory_GB} \
        -o ${id}
    cp ${id}/scaffolds.fasta ${id}.fa
    """
}

process viralverify_db_download {

    publishDir "${params.outdir}/dbs" , mode: 'copy'

    output:
    path 'nbc_hmms.hmm', emit: viralverify_db

    script:
    """
    VIRALVERIFY_DB_URL="https://ndownloader.figshare.com/files/17904323?private_link=f897d463b31a35ad7bf0"
    curl -L \${VIRALVERIFY_DB_URL} --output nbc_hmms.hmm.gz && \
        gunzip -c nbc_hmms.hmm.gz > nbc_hmms.hmm
    rm -rf nbc_hmms.hmm.gz
    """
}

process viralverify {
    tag "${id}"

    publishDir "${params.outdir}/viralverify" , mode: 'copy' ,
        pattern: "${id}/*"

    publishDir "${params.outdir}" , mode: 'copy' ,
        saveAs: { filename ->
            if (filename == "${id}/Prediction_results_fasta/${id}_plasmid.fasta") "verified_plasmids/${id}.fa"
        }

    input:
    tuple val(id), path(scaffolds)
    path(viralverify_db)

    output:
    path "${id}/*"

    script:
    """
    viralverify.py \
        -f ${scaffolds} \
        -o ${id} \
        --hmm ${viralverify_db} \
        -p \
        -t ${task.cpus}
    """
}