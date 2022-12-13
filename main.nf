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
params.query = "${projectDir}/data/query.fa" // Query genome
params.output_dir = "${projectDir}/results" // Output directory

// -----------------------------------------------------------------------------
// PROCESSES
// -----------------------------------------------------------------------------

// Align the query genome to the reference genome using Minimap2
process MINIMAP2_QUERY {

    cpus 1
    conda '/home/liam/miniconda3/envs/minimap2'
    tag "${reference_fa} x ${query_fa}"

    input:
    path(reference_fa)
    path(query_fa)

    output:
    path("query_alignment.sam"), emit: sam
    path("query_alignment.paf"), emit: paf
    
    script:
    """
    minimap2 -a -x asm5 --eqx \
        ${reference_fa} \
        ${query_fa} \
        > query_alignment.sam
    minimap2 -x asm5 --eqx \
        ${reference_fa} \
        ${query_fa} \
        > query_alignment.paf
    """

}

// Scaffold the query genome using the reference genome
process RAGTAG {

    cpus 1
    conda '/home/liam/miniconda3/envs/ragtag'
    tag "${reference_fa} x ${query_fa}"

    input:
    path(reference_fa)
    path(query_fa)

    output:
    path("ragtag_output/ragtag.scaffold.confidence.txt"), emit: confidence
    path("ragtag_output/ragtag.scaffold.fasta"), emit: fasta
    path("ragtag_output/ragtag.scaffold.stats"), emit: stats
    
    script:
    """
    ragtag.py scaffold \
        ${reference_fa} \
        ${query_fa}
    """

}

// Separate placed (scaffolded) and unplaced contigs produced by RagTag
process SEPARATE_CONTIGS {

    cpus 1
    conda '/home/liam/miniconda3/envs/biopython'
    
    input:
    path(scaffolded_query_fa)

    output:
    path("placed_seqs.fa"), emit: placed_fa
    path("unplaced_seqs.fa"), emit: unplaced_fa

    script:
    """
    python3 ${projectDir}/bin/separate-placed-unplaced-contigs.py \
        -s ${scaffolded_query_fa}
    """

}

// Align the scaffolded query genome to the reference genome
process MINIMAP2_SCAFFOLD {

    cpus 1
    conda '/home/liam/miniconda3/envs/minimap2'
    tag "${reference_fa} x ${placed_fa}"

    input:
    path(reference_fa)
    path(placed_fa)

    output:
    path("scaffold_alignment.sam"), emit: sam
    path("scaffold_alignment.paf"), emit: paf
    
    script:
    """
    minimap2 -a -x asm5 --eqx \
        ${reference_fa} \
        ${placed_fa} \
        > scaffold_alignment.sam && \
    minimap2 -x asm5 --eqx \
        ${reference_fa} \
        ${placed_fa} \
        > scaffold_alignment.paf
    """

}

// Identify structural differences between the query and reference genomes using
// SyRI
process SYRI {

    cpus 1
    conda '/home/liam/miniconda3/envs/syri'
    publishDir "${projectDir}/results/syri", mode: 'copy'
    
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
    publishDir "${projectDir}/results/plotsr", mode: 'copy'
    
    input:
    path(syri_out)
    // Need to input reference_fa and placed_fa so that plotsr can locate them
    path(reference_fa)
    path(placed_fa)

    output:
    path("plotsr_default.png"), emit: png_default
    path("plotsr_itx.png"), emit: png_itx

    script:
    """
    cat << EOF > genomes.txt
    #file	name
    ${reference_fa}	reference
    ${placed_fa}	query
    EOF

    plotsr \
        --sr ${syri_out} \
        --genomes genomes.txt \
        -o plotsr_default.png

    plotsr \
        --sr ${syri_out} \
        --genomes genomes.txt \
        -o plotsr_itx.png \
        --itx
    """

}

// Align the query genome to the reference genome using nucmer, map each
// position of each reference to its best hit in the query, and output
// alignment coordinates
process MUMMER_QUERY {

    cpus 4
    conda '/home/liam/miniconda3/envs/mummer'
    tag "${reference_fa} x ${query_fa}"    

    input:
    path(reference_fa)
    path(query_fa)

    output:
    path("query_alignment.delta"), emit: delta   
    path("query_alignment.delta.filter"), emit: filter
    path("query_alignment.delta.filter.coords"), emit: coords

    script:
    """
    nucmer -p query_alignment ${reference_fa} ${query_fa} && \
    delta-filter -r query_alignment.delta > query_alignment.delta.filter && \
    show-coords -c query_alignment.delta.filter > query_alignment.delta.filter.coords
    """

}

// Generate interactive dotplots of the genome alignments
process DOTPLOTLY {

    cpus 4
    conda '/home/liam/miniconda3/envs/dotplotly'
    publishDir "${projectDir}/results/dotplotly", mode: 'copy'

    input:
    path(nucmer_query_coords)
    path(minimap2_query_paf)
    path(minimap2_scaffold_paf)

    output:
    path("nucmer_query_dotplot.png"), emit: nucmer_query_png
    path("minimap2_query_dotplot.png"), emit: minimap2_query_png
    path("minimap2_scaffold_dotplot.png"), emit: minimap2_scaffold_png
    path("nucmer_query_dotplot.html"), emit: nucmer_query_html
    path("minimap2_query_dotplot.html"), emit: minimap2_query_html
    path("minimap2_scaffold_dotplot.html"), emit: minimap2_scaffold_html 

    script:
    """
    Rscript ${projectDir}/bin/dotPlotly/mummerCoordsDotPlotly.R -s -t -l \
        -m 0 -q 0 \
        -i ${nucmer_query_coords} \
        -o nucmer_query_dotplot && \
    Rscript ${projectDir}/bin/dotPlotly/pafCoordsDotPlotly.R -s -t -l \
        -m 0 -q 0 \
        -i ${minimap2_query_paf} \
        -o minimap2_query_dotplot && \
    Rscript ${projectDir}/bin/dotPlotly/pafCoordsDotPlotly.R -s -t -l \
        -m 0 -q 0 \
            -i ${minimap2_scaffold_paf} \
            -o minimap2_scaffold_dotplot
    """

}

// -----------------------------------------------------------------------------
// CHANNELS
// -----------------------------------------------------------------------------

reference_ch = channel.fromPath("${params.reference}", checkIfExists: true)
query_ch = channel.fromPath("${params.query}", checkIfExists: true)

// -----------------------------------------------------------------------------
// WORFLOW
// -----------------------------------------------------------------------------

workflow {

    // reference_ch.view()
    // query_ch.view()

    minimap2_query = MINIMAP2_QUERY(reference_ch, query_ch)

    mummer_query = MUMMER_QUERY(reference_ch, query_ch)

    ragtag = RAGTAG(reference_ch, query_ch)

    separate_contigs = SEPARATE_CONTIGS(ragtag.fasta)

    minimap2_scaffold = MINIMAP2_SCAFFOLD(reference_ch, separate_contigs.placed_fa)

    syri = SYRI(minimap2_scaffold.sam, reference_ch, separate_contigs.placed_fa)

    plotsr = PLOTSR(syri.out, reference_ch, separate_contigs.placed_fa)

    dotplotly = DOTPLOTLY(mummer_query.coords, minimap2_query.paf, minimap2_scaffold.paf)

}