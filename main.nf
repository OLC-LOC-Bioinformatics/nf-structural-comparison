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
process MINIMAP2_DRAFT {

    cpus 1
    conda '/home/liam/miniconda3/envs/minimap2'
    tag "${reference_fa} x ${draft_fa}"

    input:
    path(reference_fa)
    path(draft_fa)

    output:
    path("draft_alignment.sam"), emit: sam
    
    script:
    """
    minimap2 -a -x asm5 --eqx \
        ${reference_fa} \
        ${draft_fa} \
        > draft_alignment.sam
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

// Align the scaffolded draft genome to the reference genome
process MINIMAP2_SCAFFOLD {

    cpus 1
    conda '/home/liam/miniconda3/envs/minimap2'
    tag "${reference_fa} x ${placed_fa}"

    input:
    path(reference_fa)
    path(placed_fa)

    output:
    path("scaffold_alignment.sam"), emit: sam
    
    script:
    """
    minimap2 -a -x asm5 --eqx \
        ${reference_fa} \
        ${placed_fa} \
        > scaffold_alignment.sam
    """

}

// Identify structural differences between the draft and reference genomes using
// SyRI
process SYRI {

    cpus 1
    conda '/home/liam/miniconda3/envs/syri'
    
    input:
    path(alignment_sam)
    path(reference_fa)
    path(placed_fa)

    output:
    path("syri.log"), emit: log
    path("syri.out"), emit: out
    path("syri.summary"), emit: summary
    path("syri.vcf"), emit: vcf

    script:
    """
    syri -F S --no-chrmatch \
        -c ${alignment_sam} \
        -r ${reference_fa} \
        -q ${placed_fa}
    """

}

// Visualize structural differences between genomes
process PLOTSR {

    cpus 1
    conda '/home/liam/miniconda3/envs/plotsr'
    
    input:
    path(syri_out)
    // Need to input reference_fa and placed_fa so that plotsr can locate them
    path(reference_fa)
    path(placed_fa)

    output:
    path("output_plot.png"), emit: png

    script:
    """
    cat << EOF > genomes.txt
    #file	name
    reference.fa	reference
    placed_seqs.fa	draft
    EOF

    plotsr \
        --sr ${syri_out} \
        --genomes genomes.txt \
        -o output_plot.png
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

    minimap2_draft = MINIMAP2_DRAFT(reference_ch, draft_ch)

    ragtag = RAGTAG(reference_ch, draft_ch)

    separate_contigs = SEPARATE_CONTIGS(ragtag.fasta)

    minimap2_scaffold = MINIMAP2_SCAFFOLD(reference_ch, separate_contigs.placed_fa)

    syri = SYRI(minimap2_scaffold.sam, reference_ch, separate_contigs.placed_fa)

    plotsr = PLOTSR(syri.out, reference_ch, separate_contigs.placed_fa)
}