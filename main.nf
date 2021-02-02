#!/usr/bin/env nextflow

Channel
    .fromFilePairs( params.reads, size: (params.single_end || params.interleaved) ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}." }
    .set { raw_reads_deinterleave_ch }

/*
 * Step 0. Deinterleave paired reads
 */
if (params.interleaved) {
    process deinterleave {      
        tag "${id}"

        input:
        tuple val(id), path(reads) from raw_reads_deinterleave_ch

        output:
        tuple val(id), path("read_*.fastq.gz") into raw_reads_stats_ch, raw_reads_adapter_ch

        script:
        task_memory_GB = task.memory.toGiga()
        
        """
        reformat.sh \
            -Xmx${task_memory_GB}g \
            in=$reads \
            out1=read_1.fastq.gz \
            out2=read_2.fastq.gz \
            t=1
        """
    }
} else {
    raw_reads_deinterleave_ch.into { raw_reads_stats_ch; raw_reads_adapter_ch }
}

/*
 * Step 1. Raw reads histograms
 */
process raw_reads_stats {   
    tag "${id}"

    publishDir "${params.outdir}/raw_reads_stats/${id}" , mode: 'copy'

    input:
    tuple val(id), path(reads) from raw_reads_stats_ch

    output:
    path "*hist.txt"

    script:
    task_memory_GB = task.memory.toGiga()
    input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    """
    bbduk.sh \
        -Xmx${task_memory_GB}g \
        $input \
        bhist=bhist.txt \
        qhist=qhist.txt \
        gchist=gchist.txt \
        aqhist=aqhist.txt \
        lhist=lhist.txt \
        gcbins=auto \
        threads=${task.cpus}
    """
}
 
/*
 * Step 2.a Remove adapters
 */
if (!params.skip_adapter) {
    process adapter {
        tag "${id}"
    
        publishDir "${params.outdir}/qc/${id}" , mode: 'copy',
            pattern: "stats_adapter.txt"

        input:
        tuple val(id), path(reads) from raw_reads_adapter_ch

        output:
        tuple val(id), path("clean_adapter*.fastq.gz") into clean_reads_contaminant_ch
        path "stats_adapter.txt"

        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
        output = params.single_end ? "out=clean_adapter.fastq.gz" : "out1=clean_adapter_1.fastq.gz out2=clean_adapter_2.fastq.gz"
        """
        bbduk.sh \
            -Xmx${task_memory_GB}g \
            $input \
            $output \
            ref=adapters \
            ktrim=r \
            k=23 \
            mink=11 \
            hdist=1 \
            tpe \
            tbo \
            interleaved=f \
            stats=stats_adapter.txt \
            threads=${task.cpus}
        """
    }
} else {
    raw_reads_adapter_ch.set { clean_reads_contaminant_ch }
}

/*
 * Step 2.b Remove contaminants
 */
if (!params.skip_contaminant) {
    process contaminant {
        tag "${id}"
    
        publishDir "${params.outdir}/qc/${id}" , mode: 'copy', 
            pattern: "stats_contaminant.txt"

        input:
        tuple val(id), path(reads) from clean_reads_contaminant_ch

        output:
        tuple val(id), path("clean_contaminant*.fastq.gz") into clean_reads_quality_ch
        path "stats_contaminant.txt"

        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
        output = params.single_end ? "out=clean_contaminant.fastq.gz" : "out1=clean_contaminant_1.fastq.gz out2=clean_contaminant_2.fastq.gz"
        """
        bbduk.sh \
            -Xmx${task_memory_GB}g \
            $input \
            $output \
            ref=artifacts,phix \
            k=31 \
            hdist=1 \
            interleaved=f \
            stats=stats_contaminant.txt \
            threads=${task.cpus}
        """
    }
} else {
    clean_reads_contaminant_ch.set { clean_reads_quality_ch }
}

/*
 * Step 2.c Quality filtering/trimming and length filtering
 */
if (!params.skip_quality) {
    process quality {
        tag "${id}"

        input:
        tuple val(id), path(reads) from clean_reads_quality_ch

        output:
        tuple val(id), path("clean*.fastq.gz") into \
            clean_reads_stats_ch, clean_reads_spades_ch, clean_reads_megahit_ch
  
        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
        output = params.single_end ? "out=clean.fastq.gz" : "out1=clean_1.fastq.gz out2=clean_2.fastq.gz"
        """
        bbduk.sh \
            $input \
            $output \
            maq=10 \
            maxns=4 \
            qtrim=r \
            trimq=6 \
            mlf=0.5 \
            minlen=50 \
            interleaved=f \
            threads=${task.cpus}
            """
    }
} else {
    clean_reads_quality_ch.into {
        clean_reads_stats_ch;
        clean_reads_spades_ch;
        clean_reads_megahit_ch
        }
}

/*
 * Step 3. Clean reads histograms
 */
process clean_reads_stats {
    tag "${id}"

    publishDir "${params.outdir}/clean_reads_stats/${id}" , mode: 'copy'

    when:
    ! (params.skip_adapter && params.skip_contaminant && params.skip_quality)

    input:
    tuple val(id), path(reads) from clean_reads_stats_ch

    output:
    path "*hist.txt"

    script:
    task_memory_GB = task.memory.toGiga()
    input = params.single_end ? "in=\"$reads\"" : "in1=\"${reads[0]}\" in2=\"${reads[1]}\""
    """
    bbduk.sh \
        -Xmx${task_memory_GB}g \
        $input \
        bhist=bhist.txt \
        qhist=qhist.txt \
        gchist=gchist.txt \
        aqhist=aqhist.txt \
        lhist=lhist.txt \
        gcbins=auto \
        threads=${task.cpus}
    """
}

/*
 * Step 4.a Assembly with Spades
 */
if (!params.single_end && !params.megahit_only) {
    process spades {
        tag "${id}"
    
        publishDir "${params.outdir}" , mode: 'copy'   

        publishDir "${params.outdir}" , mode: 'copy' ,
            saveAs: { filename ->
                if (filename == "spades/scaffolds.fasta") "scaffolds/${id}.scaffolds.fa"
            }
    
        input:
        tuple val(id), path(reads) from clean_reads_spades_ch
    
        output:
        tuple val(id), path("spades/scaffolds.fasta") into scaffolds_spades_ch
        path "spades/*"
    
        script:
        task_memory_GB = task.memory.toGiga()
        param_spades_k = params.spades_k == 'default' ? "" :  "-k ${params.spades_k}"
        """
        spades.py \
            --isolate \
            -1 ${reads[0]} \
            -2 ${reads[1]} \
            ${param_spades_k} \
            --threads ${task.cpus} \
            --memory ${task_memory_GB} \
            -o spades
        """
    }
} else {
    scaffolds_spades_ch = Channel.empty()
}

/*
 * Step 4.b Assembly with Megahit
 */
if (params.single_end || params.megahit_only) {
    process megahit {
        tag "${id}"

        publishDir "${params.outdir}" , mode: 'copy'

        publishDir "${params.outdir}" , mode: 'copy' ,
            saveAs: { filename -> 
                if (filename == "megahit/final.contigs.fa") "scaffolds/${id}.scaffolds.fa"
            }

        input:
        tuple val(id), path(reads) from clean_reads_megahit_ch

        output:
        tuple val(id), path("megahit/final.contigs.fa") into scaffolds_megahit_ch
        path "megahit/*"

        script:
        task_memory_GB = task.memory.toGiga()
        input = params.single_end ? "-r \"$reads\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
        param_megahit_k = params.megahit_k == 'default' ? "" :  "--k-list ${params.megahit_k}"
        """
        megahit \
            $input \
            -t ${task.cpus} \
            ${param_megahit_k} \
            --min-count 3 \
            --memory $task_memory_GB \
            -o megahit
        """
    }
} else {
    scaffolds_megahit_ch = Channel.empty()
}

scaffolds_spades_ch
    .mix(scaffolds_megahit_ch)
    .into { scaffolds_stats_ch }

/*
 * Step 5. Scaffold statistics
 */
process statswrapper {      
        tag "${id}"

        publishDir "${params.outdir}" , mode: 'copy'

        input:
        tuple val(id), path(scaffolds) from scaffolds_stats_ch

        output:
        path 'stats.tsv'

        script:       
        """
        mkdir scaffolds_dir
        mv ${scaffolds} scaffolds_dir
        statswrapper.sh scaffolds_dir/* > stats.tsv
        """
}
