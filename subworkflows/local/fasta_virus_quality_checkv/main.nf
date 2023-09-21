//
// Assess virus quality with Checkv
//

include { CHECKV_DOWNLOADDATABASE   } from '../../../modules/nf-core/checkv/downloaddatabase/main'
include { GUNZIP                    } from '../../../modules/nf-core/gunzip/main'
include { CHECKV_ENDTOEND           } from '../../../modules/nf-core/checkv/endtoend/main'

workflow FASTA_VIRUS_QUALITY_CHECKV {
    take:
    virus_fasta_gz  // [ [ meta ], fasta.gz ]       , assemblies/genomes (mandatory)
    checkv_db       // [ [ meta ], checkv_db_dir ]  , CheckV database directory (optional)

    main:
    ch_versions = Channel.empty()

    // if genomad_db exists, skip CHECKV_DOWNLOADDATABASE
    if ( checkv_db ){
        ch_checkv_db = checkv_db
    } else {
        //
        // MODULE: download CheckV database
        //
        ch_checkv_db_dir = CHECKV_DOWNLOADDATABASE( ).checkv_db
        ch_checkv_db = ch_checkv_db_dir.map{[ [ id:'checkv_db' ], it ]}
        ch_versions = ch_versions.mix(CHECKV_DOWNLOADDATABASE.out.versions.first())
    }

    //
    // MODULE: Gunzip fasta for checkV
    //
    ch_viruses_fasta = GUNZIP ( virus_fasta_gz ).gunzip
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())

    //
    // MODULE: Assess virus quality
    //
    CHECKV_ENDTOEND ( virus_fasta_gz, ch_checkv_db )
    ch_quality_summary_tsv = CHECKV_ENDTOEND.out.quality_summary
    ch_viruses_fasta_gz = CHECKV_ENDTOEND.out.viruses
    ch_proviruses_fasta_gz = CHECKV_ENDTOEND.out.proviruses
    ch_versions = ch_versions.mix(CHECKV_ENDTOEND.out.versions.first())

    emit:
    viruses_fasta_gz    = ch_viruses_fasta_gz       // [ [ meta ], viruses.fasta.gz ]       , FASTA file containing viruses
    proviruses_fasta_gz = ch_proviruses_fasta_gz    // [ [ meta ], proviruses.fasta.gz ]    , FASTA file containing proviruses
    quality_summary_tsv = ch_quality_summary_tsv    // [ [ meta ], quality_summary.tsv ]    , TSV file containing quality data
    versions            = ch_versions               // [ versions.yml ]
}
