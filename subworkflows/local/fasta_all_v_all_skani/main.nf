//
// Compare sequences by performing all-v-all skani
//
include { SKANI_SKETCH      } from '../../../modules/local/skani/sketch/main'
include { SKANI_TRIANGLE    } from '../../../modules/local/skani/triangle/main'

workflow FASTA_ALL_V_ALL_SKANI {

    take:
    fasta_gz    // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Make skani sketch
    //
    ch_sketch = SKANI_SKETCH ( fasta_gz ).db
    ch_versions = ch_versions.mix ( SKANI_SKETCH.out.versions )

    //
    // MODULE: Perform BLASTN all-v-all alignment
    //
    ch_skani_tsv    = SKANI_TRIANGLE ( ch_sketch , ch_sketch ).txt
    ch_versions     = ch_versions = ch_versions.mix ( SKANI_TRIANGLE.out.versions )

    emit:
    skani_tsv   = ch_skani_tsv  // [ [ meta ], ani.tsv ]    , TSV file containing skani all-v-all ANI and AF values
    ch_versions = ch_versions   // [ versions.yml ]

}
