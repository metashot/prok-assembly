# prok-assembly

metashot/prok-assembly is a workflow for the assembly of Illumina sequences from
prokaryotic isolates and other small genomes.

- [MetaShot Home](https://metashot.github.io/)

## Main features

- Input: single-end, paired-end (also interleaved) Illumina sequences (gzip and
  bzip2 compressed FASTQ also supported);
- Histogram text files (for each input sample) of base frequency, quality
  scores, GC content, average quality and length are generated from input reads
  and clean reads using
  [bbduk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/);
- Adapter trimming, contaminant filtering and a conservative quality
  filtering/trimming and length filtering using bbduk;
- Optionally, assemble plasmids with plasmidSpades and verify them using
  ViralVerify.

## Quick start

1. Install Docker (or Singulariry) and Nextflow (see
   [Dependences](https://metashot.github.io/#dependencies));
1. Start running the analysis:
   
  ```bash
  nextflow run metashot/prok-assembly \
    --reads '*_R{1,2}.fastq.gz' \
    --outdir results
  ```

## Parameters
See the file [`nextflow.config`](nextflow.config) for the complete list of
parameters.

## Output
The files and directories listed below will be created in the `results` directory
after the pipeline has finished.

### Main outputs
- `scaffolds`: scaffolds for each input samples;
- `stats_scaffolds.tsv`: scaffold statistics;
- `verified_plasmids`: verified plasmids (if `--run_plasmidspades` is set).

### Secondary outputs

- `raw_reads_stats`: base frequency, quality scores, gc content, average
  quality and length for each input sample;
- `clean_reads_stats`: same as above, but for the reads after the quality
  control;
- `clean_reads`: clean reads (if `--save_clean` is set);
- `qc`: statistics about the adapter trimming and the contaminant filtering;
- `spades`: SPAdes output for each input sample (if `--save_spades_output` is
  set);
- `scaffolds_plasmids`: candidate plasmids (if `--run_plasmidspades` is set);
- `viralverify`: viralVerify output (if `--run_plasmidspades` is set);


## System requirements
Please refer to [System
requirements](https://metashot.github.io/#system-requirements) for the complete
list of system requirements options.