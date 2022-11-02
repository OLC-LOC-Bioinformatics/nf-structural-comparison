#!/usr/bin/env nextflow

/*

    structural-comparison
    ----------
    A DSL2 Nextflow workflow for the structural comparison of bacterial genomes.

*/ 

// -----------------------------------------------------------------------------
// GLOBAL PARAMETERS
// -----------------------------------------------------------------------------

params.reference = "${projectDir}/data/reference.fa" // Reference genome
params.draft = "${projectDir}/data/draft.fa" // Draft genome
params.output_dir = "${projectDir}/results" // Output directory

// -----------------------------------------------------------------------------
// PROCESSES
// -----------------------------------------------------------------------------

// Align the draft genome to the reference genome
process MINIMAP2 {

    cpus 1
    conda '/home/liam/miniconda3/envs/minimap2'
    tag "${reference_fa} x ${draft_fa}"

    input:
    path(reference_fa)
    path(draft_fa)

    output:
    path("alignment.sam"), emit: sam
    
    script:
    """
    minimap2 -a -x asm5 --eqx \
        ${reference_fa} \
        ${draft_fa} \
        > alignment.sam
    """

}

// Scaffold the draft genome using the reference genome
process RAGTAG {

    cpus 1
    conda '/home/liam/miniconda3/envs/ragtag'
    tag "${reference_fa} x ${draft_fa}"

    input:
    path(reference_fa)
    path(draft_fa)

    output:
    path("ragtag_output/ragtag.scaffold.confidence.txt"), emit: confidence
    path("ragtag_output/ragtag.scaffold.fasta"), emit: fasta
    path("ragtag_output/ragtag.scaffold.stats"), emit: stats
    
    script:
    """
    ragtag.py scaffold \
        ${reference_fa} \
        ${draft_fa}
    """

}

// Separate placed (scaffolded) and unplaced contigs produced by RagTag
process SEPARATE_CONTIGS {

    cpus 1
    conda '/home/liam/miniconda3/envs/biopython'
    
    input:
    path(scaffolded_draft_fa)

    output:
    path("placed_seqs.fa"), emit: placed_fa
    path("unplaced_seqs.fa"), emit: unplaced_fa

    script:
    """
    python3 ${projectDir}/bin/separate-placed-unplaced-contigs.py \
        -s ${scaffolded_draft_fa}
    """

}


// -----------------------------------------------------------------------------
// CHANNELS
// -----------------------------------------------------------------------------

reference_ch = channel.fromPath("${params.reference}", checkIfExists: true)
draft_ch = channel.fromPath("${params.draft}", checkIfExists: true)

// -----------------------------------------------------------------------------
// WORFLOW
// -----------------------------------------------------------------------------

workflow {

    // reference_ch.view()
    // draft_ch.view()

    align_genomes = MINIMAP2(reference_ch, draft_ch)
    // align_genomes.view()

    scaffold_draft = RAGTAG(reference_ch, draft_ch)
    // scaffold_draft.fasta.view()

    separate_contigs = SEPARATE_CONTIGS(scaffold_draft.fasta)

}