//
// Identify sequences contained in readset
//

include { MASH_SKETCH as MASH_SKETCH_ASSEMBLIES } from '../../../modules/nf-core/mash/sketch/main'      // TODO: Update nf-core module to remove optional -r argument
include { MASH_SKETCH as MASH_SKETCH_REFERENCES } from '../../../modules/nf-core/mash/sketch/main'      // TODO: Update nf-core module to remove optional -r argument
include { MASH_PASTE                            } from '../../../modules/local/mash/paste/main'         // TODO: Add module to nf-core
include { MASH_SCREEN                           } from '../../../modules/nf-core/mash/screen/main'      // TODO: Update nf-core module to add meta to sketch

workflow FASTQ_FASTA_REFERENCE_CONTAINMENT_MASH {
    take:
    fastq_gz                // [ [ meta.id ] , [ 1.fastq.gz, 2.fastq.gz ]  , reads (optional)
    assembly_fasta_gz       // [ [ meta.id ] , fasta.gz ]                  , assemblies (optional)
    reference_fasta_gz      // [ [ meta.id ] , fasta.gz ]                  , reference sequences (mandatory)
    reference_sketch_msh    // [ [ meta.id ] , fasta.msh ]                 , reference_sequences_sketch (optional)

    main:
    ch_versions = Channel.empty()

    // if provided, use reference sketch. If not, create one
    if ( reference_sketch_msh ) {
        ch_reference_sketch_msh = reference_sketch_msh
    } else {
        //
        // MODULE: Create sketch of reference sequences
        //
        ch_reference_sketch_msh = MASH_SKETCH_REFERENCES ( reference_fasta_gz ).mash
        ch_versions = ch_versions.mix(MASH_SKETCH_REFERENCES.out.versions.first())
    }

    // identify which assemblies have an associated fastq file
    ch_fasta_gz_w_fastq_gz = fastq_gz
        .join( assembly_fasta_gz, by:0 )
        .map {
            meta, fastq_gz, fasta_gz ->
            [ meta, fasta_gz]
        }

    //
    // MODULE: Create a sketch of assemblies that have associated fastq files
    //
    ch_assembly_sketch_msh = MASH_SKETCH_ASSEMBLIES ( ch_fasta_gz_w_fastq_gz ).mash
    ch_versions = ch_versions.mix(MASH_SKETCH_ASSEMBLIES.out.versions.first())

    //
    // MODULE: Combine assembly and reference sketches
    //
    ch_combined_sketch_msh = MASH_PASTE ( ch_assembly_sketch_msh, ch_reference_sketch_msh.collect() ).msh
    ch_versions = ch_versions.mix(MASH_PASTE.out.versions.first())

    // split fastq files based on whether or not they have a paired assembly file
    ch_fastq_gz_branch = fastq_gz
        .join( assembly_fasta_gz, by:0, remainder: true )
        .branch {
            meta, fastq_gz, fasta_gz ->
                fastq_fasta_pair: fastq_gz.size() > 0 && fasta_gz
                fastq_only: fastq_gz.size() > 0 && !fasta_gz
        }
    ch_fastq_gz_fasta_pair = ch_fastq_gz_branch.fastq_fasta_pair
        .map {
            meta, fastq_gz, fasta_gz ->
            [ meta, fastq_gz ]
        }
    ch_fastq_gz_only = ch_fastq_gz_branch.fastq_only
        .map {
            meta, fastq_gz, fasta_gz ->
            [ meta, fastq_gz ]
        }

    ch_fastq_gz_branch.fastq_fasta_pair.view()

    // join reads and combined sketch by meta.id if assemblies are included
    ch_screen_input_for_paired = ch_fastq_gz_fasta_pair.join( ch_combined_sketch_msh, by:0 )

    // join reads and reference sketch if no assemblies are included
    ch_screen_input_for_fastq_only = ch_fastq_gz_only.combine( ch_reference_sketch_msh.map { it[1] } )

    // combine the two options for input into mash screen
    ch_mash_screen_input = ch_screen_input_for_paired.mix( ch_screen_input_for_fastq_only )

    //
    // MODULE: Identify contained genomes
    //
    ch_mash_screen_tsv = MASH_SCREEN ( ch_mash_screen_input ).screen
    ch_versions = ch_versions.mix(MASH_SCREEN.out.versions.first())

    emit:
    mash_screen_results = ch_mash_screen_tsv    // [ [ meta.id ] , fasta.gz ]   , concatenated assemblies and contained references
    versions            = ch_versions           // [ versions.yml ]
}
