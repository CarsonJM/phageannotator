//
// Assess virus quality with Checkv
//
include { CHECKV_DOWNLOADDATABASE   } from '../../../modules/nf-core/checkv/downloaddatabase/main'
include { CHECKV_ENDTOEND           } from '../../../modules/nf-core/checkv/endtoend/main'
include { QUALITYFILTERVIRUSES      } from '../../../modules/local/qualityfilterviruses/main'

workflow FASTA_VIRUS_QUALITY_CHECKV {
    take:
    virus_fasta_gz  // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)
    checkv_db       // [ checkv_db ]            , CheckV database directory (optional)

    main:
    ch_versions = Channel.empty()

    // if checkv DB exists, skip CHECKV_DOWNLOADDATABASE
    if ( checkv_db ){
        ch_checkv_db = checkv_db
    } else {
        //
        // MODULE: download standard CheckV database
        //
        ch_checkv_db    = CHECKV_DOWNLOADDATABASE( ).checkv_db
        ch_versions     = ch_versions.mix( CHECKV_DOWNLOADDATABASE.out.versions )
    }

    //
    // MODULE: Run CheckV end to end
    //
    ch_quality_summary_tsv  = CHECKV_ENDTOEND ( virus_fasta_gz, ch_checkv_db.collect() )
    ch_viruses_fna_gz       = CHECKV_ENDTOEND.out.viruses
    ch_proviruses_fna_gz    = CHECKV_ENDTOEND.out.proviruses
    ch_versions             = ch_versions.mix( CHECKV_ENDTOEND.out.versions )

    // create channel for input into QUALITY_FILTER_VIRUSES
    ch_quality_filter_viruses_input = ch_viruses_fna_gz
        .join( ch_proviruses_fna_gz )
        .join( ch_quality_summary_tsv )
        .multiMap { it ->
            viruses: [ it[0], it[1] ]
            proviruses: [ it[0], it[2] ]
            quality_summary: [ it[0], it[3] ]
        }

    //
    // MODULE: Quality filter viruses
    //
    ch_filtered_viruses_fna_gz = QUALITYFILTERVIRUSES (
        ch_quality_filter_viruses_input.viruses,
        ch_quality_filter_viruses_input.proviruses,
        ch_quality_filter_viruses_input.quality_summary
    ).filtered_viruses
    ch_versions = ch_versions.mix( QUALITYFILTERVIRUSES.out.versions )

    emit:
    filtered_viruses_fna_gz     = ch_quality_filtered_viruses_fna_gz    // [ [ meta ], viruses.fna.gz ]       , FASTA file containing quality-filtered viruses
    quality_summary_tsv         = ch_quality_summary_tsv                // [ [ meta ], quality_summary.tsv ]  , TSV file containing quality data
    versions                    = ch_versions                           // [ versions.yml ]
}
