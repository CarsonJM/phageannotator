//
// Compare sequences by performing an all-v-all BLAST
//
include { BLAST_MAKEBLASTDB         } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN              } from '../../../modules/nf-core/blast/blastn/main'
include { CHECKVANICLUSTER_ANICALC  } from '../../../modules/local/checkvanicluster/anicalc/main'

workflow FASTA_ALL_V_ALL_BLAST {

    take:
    fasta_gz    // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Make BLASTN database
    //
    ch_blast_db = BLAST_MAKEBLASTDB ( fasta_gz ).db
    ch_versions = ch_versions.mix( BLAST_MAKEBLASTDB.out.versions )

    //
    // MODULE: Perform BLASTN all-v-all alignment
    //
    ch_blast_txt = BLAST_BLASTN ( fasta_gz , ch_blast_db ).txt
    ch_versions = ch_versions = ch_versions.mix( BLAST_MAKEBLASTDB.out.versions )

    //
    // MODULE: Calculate average nucleotide identity (ANI) and alignment fraction (AF) based on BLAST
    //
    ch_blast_ani_tsv    = CHECKVANICLUSTER_ANICALC ( ch_blast_txt ).ani
    ch_versions         = ch_versions.mix( CHECKVANICLUSTER_ANICALC.out.versions )

    emit:
    blast_ani_tsv   = ch_blast_ani_tsv  // [ [ meta ], ani.tsv ]   , TSV file containing BLAST all-v-all ANI and AF values
    ch_versions     = ch_versions       // [ versions.yml ]

}
