/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULES: Local modules
//
include { SEQKIT_SEQ                                } from '../../modules/local/seqkit/seq/main'                                    // TODO: Add to nf-core
include { EXTRACTVIRALASSEMBLIES                    } from '../../modules/local/extractviralassemblies/main'
include { QUALITYFILTERVIRUSES                      } from '../../modules/local/qualityfilterviruses/main'
include { SKANI_TRIANGLE                            } from '../../modules/local/skani/triangle/main'
include { COVERM_CONTIG                             } from '../../modules/local/coverm/contig/main'                                 // TODO: Add to nf-core
include { INSTRAIN_STB                              } from '../../modules/local/instrain/stb/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { methodsDescriptionText                    } from '../../subworkflows/local/utils_nfcore_phageannotator_pipeline'
include { FASTQ_VIRUS_ENRICHMENT_VIROMEQC           } from '../../subworkflows/local/fastq_virus_enrichment_viromeqc/main'
include { FASTA_VIRUS_CLASSIFICATION_GENOMAD        } from '../../subworkflows/local/fasta_virus_classification_genomad/main'       // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTQ_FASTA_CONTIG_EXTENSION_COBRA        } from '../../subworkflows/local/fastq_fasta_contig_extension_cobra/main'       // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_VIRUS_QUALITY_CHECKV                } from '../../subworkflows/local/fasta_virus_quality_checkv/main'               // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_ALL_V_ALL_BLAST                     } from '../../subworkflows/local/fasta_all_v_all_blast/main'
include { FASTA_PHAGE_HOST_IPHOP                    } from '../../subworkflows/local/fasta_phage_host_iphop/main'                   // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTA_PHAGE_FUNCTION_PHAROKKA             } from '../../subworkflows/local/fasta_phage_function_pharokka/main'
include { FASTA_MICRODIVERSITY_INSTRAIN             } from '../../subworkflows/local/fasta_microdiversity_instrain/main'            // TODO: Add to nf-core; Add nf-tests to nf-core modules
include { FASTQ_FASTA_PROVIRUS_ACTIVITY_PROPAGATE   } from '../../subworkflows/local/fastq_fasta_provirus_activity_propagate/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// PLUGINS
//
include { paramsSummaryMap       } from 'plugin/nf-validation'

//
// MODULES: Installed directly from nf-core/modules
//
include { BACPHLIP                              } from '../../modules/nf-core/bacphlip/main'
include { CAT_CAT as CAT_VIRUSES                } from '../../modules/nf-core/cat/cat/main'
include { FASTQC                                } from '../../modules/nf-core/fastqc/main'
include { GENOMAD_ENDTOEND as GENOMAD_TAXONOMY  } from '../../modules/nf-core/genomad/endtoend/main'
include { MULTIQC                               } from '../../modules/nf-core/multiqc/main'
include { NUCMER                                } from '../../modules/nf-core/nucmer/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHAGEANNOTATOR {

    take:
    fastq_gz    // [ [ meta ], reads.fastq.gz ]     , reads (mandatory)
    fasta_gz    // [ [ meta ], assembly.fasta.gz ]  , assemblies/genomes (mandatory)

    main:
    ch_multiqc_files    = Channel.empty()
    ch_versions         = Channel.empty()

    // make a channel for each item that contains a fastq file
    ch_fastq_gz = fastq_gz
        .branch {
            meta, fastq_gz ->
                fastq_included: fastq_gz[0].size() > 0
                fastq_missing: fastq_gz[0].size() == 0
        }

    //
    // MODULE: Run FastQC
    //
    FASTQC ( ch_fastq_gz.fastq_included )
    ch_multiqc_files    = ch_multiqc_files.mix ( FASTQC.out.zip.collect{it[1]} )
    ch_versions         = ch_versions.mix ( FASTQC.out.versions  )


    /*----------------------------------------------------------------------------
        Estimate viral enrichment in reads
    ------------------------------------------------------------------------------*/
    if ( params.run_viromeqc ) {
        //
        // SUBWORKFLOW: Download and run ViromeQC
        //
        ch_virus_enrichment_tsv = FASTQ_VIRUS_ENRICHMENT_VIROMEQC ( ch_fastq_gz.fastq_included ).enrichment_tsv
        ch_versions             = ch_versions.mix ( FASTQ_VIRUS_ENRICHMENT_VIROMEQC.out.versions )
    } else {
        ch_virus_enrichment_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Remove low-length assemblies
    ------------------------------------------------------------------------------*/
    if ( params.run_seqkit_seq ) {
        //
        // MODULE: Filter assemblies by length
        //
        ch_filtered_input_fasta_gz  = SEQKIT_SEQ ( fasta_gz ).fastx
        ch_versions                 = ch_versions.mix ( SEQKIT_SEQ.out.versions )
    } else {
        ch_filtered_input_fasta_gz  = fasta_gz
    }


    /*----------------------------------------------------------------------------
        De novo virus classification
    ------------------------------------------------------------------------------*/
    if ( params.run_genomad || params.run_genomad_taxonomy ) {
        // create channel from params.genomad_db
        if ( !params.genomad_db ){
            ch_genomad_db = null
        } else {
            ch_genomad_db = Channel.value( file( params.genomad_db, checkIfExists:true ) )
        }

        if ( params.run_genomad ) {
            //
            // SUBWORKFLOW: Download and run geNomad
            //
            ch_viruses_fna_gz       = FASTA_VIRUS_CLASSIFICATION_GENOMAD ( ch_filtered_input_fasta_gz, ch_genomad_db ).viruses_fna_gz
            ch_genomad_db_dir       = FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.genomad_db
            ch_versions             = ch_versions.mix ( FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.versions )
            ch_virus_summaries_tsv  = FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.virus_summaries_tsv
        } else {
            //
            // SUBWORKFLOW: Download genomad database
            //
            ch_genomad_db_dir       = FASTA_VIRUS_CLASSIFICATION_GENOMAD ( Channel.empty(), ch_genomad_db ).genomad_db
            ch_versions             = ch_versions.mix ( FASTA_VIRUS_CLASSIFICATION_GENOMAD.out.versions )
            ch_viruses_fna_gz       = ch_filtered_input_fasta_gz
            ch_virus_summaries_tsv  = Channel.empty()
        }
    } else {
        ch_viruses_fna_gz       = ch_filtered_input_fasta_gz
        ch_virus_summaries_tsv  = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Extend viral contigs
    ------------------------------------------------------------------------------*/
    if ( params.run_cobra ) {
        //
        // SUBWORKFLOW: Extend assembled contigs
        //
        ch_extended_viruses_fasta_gz    = FASTQ_FASTA_CONTIG_EXTENSION_COBRA (
            ch_fastq_gz.fastq_included,
            fasta_gz,
            ch_virus_summaries_tsv,
            params.cobra_assembler,
            params.cobra_mink,
            params.cobra_maxk
        ).extended_fasta
        ch_virus_extension_summary_tsv  = FASTQ_FASTA_CONTIG_EXTENSION_COBRA.out.cobra_summary_tsv
        ch_versions                     = ch_versions.mix( FASTQ_FASTA_CONTIG_EXTENSION_COBRA.out.versions )
    } else {
        ch_extended_viruses_fasta_gz    = ch_viruses_fna_gz
        ch_virus_extension_summary_tsv  = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Assess virus quality and filter
    ------------------------------------------------------------------------------*/
    if ( params.run_checkv ) {
        // create channel from params.checkv_db
        if ( !params.checkv_db ){
            ch_checkv_db = null
        } else {
            ch_checkv_db = Channel.value( file( params.checkv_db, checkIfExists:true ) )
        }

        //
        // SUBWORKFLOW: Assess virus quality
        //
        ch_quality_summary_tsv      = FASTA_VIRUS_QUALITY_CHECKV ( ch_extended_viruses_fasta_gz, ch_checkv_db ).quality_summary_tsv
        ch_filtered_viruses_fna_gz  = FASTA_VIRUS_QUALITY_CHECKV.out.filtered_viruses_fna_gz
        ch_versions                 = ch_versions.mix( FASTA_VIRUS_QUALITY_CHECKV.out.versions )
    } else {
        ch_filtered_viruses_fna_gz  = ch_extended_viruses_fasta_gz
        ch_quality_summary_tsv      = Channel.empty()
    }

    /*----------------------------------------------------------------------------
        Assign viral taxonomy
    ------------------------------------------------------------------------------*/
    if ( params.run_genomad_taxonomy ) {
        //
        // SUBWORKFLOW: Assign taxonomy using ICTV taxa specific marker genes
        //
        ch_marker_taxonomy_tsv  = GENOMAD_TAXONOMY ( ch_filtered_viruses_fna_gz, ch_genomad_db_dir ).taxonomy
        ch_versions             = ch_versions.mix( GENOMAD_TAXONOMY.out.versions )
    } else {
        ch_marker_taxonomy_tsv  = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Predict phage hosts
    ------------------------------------------------------------------------------*/
    if ( params.run_iphop ){
        // create channel from params.checkv_db
        if ( !params.iphop_db ){
            ch_iphop_db = null
        } else {
            ch_iphop_db = file( params.iphop_db, checkIfExists:true )
        }

        //
        // SUBWORKFLOW: Download database and predict phage hosts
        //
        ch_host_predictions_tsv = FASTA_PHAGE_HOST_IPHOP ( ch_filtered_viruses_fna_gz, ch_iphop_db ).host_predictions_tsv
        ch_versions             = ch_versions.mix( FASTA_PHAGE_HOST_IPHOP.out.versions )
    } else {
        ch_host_predictions_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Phage functional annotation
    ------------------------------------------------------------------------------*/
    if ( params.run_pharokka ) {
        // create channel from params.pharokka_db
        if ( !params.pharokka_db ){
            ch_pharokka_db = null
        } else {
            ch_pharokka_db = Channel.value( file( params.pharokka_db, checkIfExists:true ) )
        }

        //
        // SUBWORKFLOW: Functionally annotate phage sequences
        //
        ch_pharokka_gbk_gz      = FASTA_PHAGE_FUNCTION_PHAROKKA ( ch_anicluster_reps_fasta, ch_pharokka_db ).pharokka_gbk_gz
        ch_pharokka_output_tsv  = FASTA_PHAGE_FUNCTION_PHAROKKA.out.pharokka_final_output_tsv
        ch_versions             = ch_versions.mix( FASTA_PHAGE_FUNCTION_PHAROKKA.out.versions )
    } else {
        ch_pharokka_gbk_gz      = Channel.empty()
        ch_pharokka_output_tsv  = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Predict virus lifestyle
    ------------------------------------------------------------------------------*/
    if ( params.run_bacphlip ) {
        //
        // MODULE: Predict phage lifestyle with BACPHLIP
        //
        ch_bacphlip_lifestyle_tsv   = BACPHLIP ( ch_anicluster_reps_fasta ).bacphlip_results
        ch_versions                 = ch_versions.mix( BACPHLIP.out.versions )
    } else {
        ch_bacphlip_lifestyle_tsv   = Channel.empty()
    }

    // TODO: Add integration status check
    // TODO: Add PHROG integrase identification


    /*----------------------------------------------------------------------------
        Predict if proviruses are active
    ------------------------------------------------------------------------------*/
    // TODO: Update python code with Adam's recommendations
    if ( params.run_propagate ){
        ch_provirus_activity_tsv    = FASTQ_FASTA_PROVIRUS_ACTIVITY_PROPAGATE (
            ch_fastq_gz.fastq_included,
            fasta_gz,
            ch_virus_summaries_tsv,
            ch_quality_summary_tsv,
            ch_clusters_tsv,
            params.propagate_min_ani,
            params.propagate_min_qcov,
            params.propagate_min_tcov
            ).propagate_results_tsv
        ch_versions                 = ch_versions.mix ( FASTQ_FASTA_PROVIRUS_ACTIVITY_PROPAGATE.out.versions )
    } else {
        ch_provirus_activity_tsv    = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Calculate all-v-all ANI
    ------------------------------------------------------------------------------*/
    if ( params.run_blast_ani || params.run_skani_ani || params.run_nucmer_ani ) {
        // create a channel for combining filtered viruses (sorted so output is the same for tests)
        ch_cat_viruses_input = ch_filtered_viruses_fna_gz
                                .map { [ [ id:'all_samples' ], it[1] ] }
                                .groupTuple( sort: 'deep' )

        //
        // MODULE: Concatenate all quality filtered viruses into one file
        //
        ch_filtered_viruses_combined_fna_gz = CAT_VIRUSES ( ch_cat_viruses_input ).file_out
        ch_versions                         = ch_versions.mix( CAT_VIRUSES.out.versions )

        if ( params.run_blast_ani ){
            //
            // SUBWORKFLOW: Calculate BLAST all-v-all based ANI
            //
            ch_blast_ani_tsv    = FASTA_ALL_V_ALL_BLAST ( ch_filtered_viruses_combined_fna_gz ).blast_ani_tsv
            ch_versions         = ch_versions.mix( FASTA_ALL_V_ALL_BLAST.out.versions )
        }
        if ( params.run_skani_ani ) {
            //
            // SUBWORKFLOW: Calculate skani all-v-all based ANI
            //
            ch_skani_ani_tsv    = SKANI_TRIANGLE ( ch_filtered_viruses_combined_fna_gz ).ani_tsv
            ch_versions         = ch_versions.mix( SKANI_TRIANGLE.out.versions )
        }
        if ( params.run_nucmer_ani ) {
            //
            // SUBWORKFLOW: Calculate skani all-v-all based ANI
            //
            ch_nucmer_ani_tsv   = NUCMER ( ch_filtered_viruses_combined_fna_gz ).coords
            ch_versions         = ch_versions.mix( NUCMER.out.versions )
        }
    } else {
        ch_blast_ani_tsv        = Channel.empty()
        ch_skani_ani_tsv        = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Calculate all-v-all AAI
    ------------------------------------------------------------------------------*/
    if ( params.run_diamond_aai || params.run_mmseqs2_aai ) {
        if ( params.run_diamond_aai ){
            //
            // SUBWORKFLOW: Calculate DIAMOND all-v-all based ANI
            //
        } else if ( params.run_mmseqs2_aai ) {
            //
            // SUBWORKFLOW: Calculate MMSeqs2 all-v-all based ANI
            //
        }
    } else {
        ch_diamond_aai_tsv      = Channel.empty()
        ch_mmseqs2_aai_tsv      = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Align reads to viruses
    ------------------------------------------------------------------------------*/
    if ( params.run_coverm || params.run_instrain ) {
        //
        // MODULE: Calculate abundance metrics with CoverM
        //
        ch_coverm_tsv   = COVERM_CONTIG ( ch_fastq, ch_dereplicated_viruses ).alignment_results
        ch_coverm_bam   = COVERM_CONTIG.out.bam
        ch_versions     = ch_versions.mix( COVERM_CONTIG.out.versions )
    } else {
        ch_coverm_tsv   = Channel.empty()
        if ( params.run_instrain ) {
            error "[nf-core/phageannotator] ERROR: skip_read_alignment = true but skip_instrain = false; read alignment must take place for inStrain to run"
        }
    }


    /*----------------------------------------------------------------------------
        Analyze phage microdiversity
    ------------------------------------------------------------------------------*/
    // TODO: Add option to run instrain/compare within groups rather than across all samples
    // TODO: Add prodigal-gv to predict genes
    if ( params.run_instrain ) {
        //
        // SUBWORKFLOW: Assess virus microdiversity within and across samples
        //
    } else {
        ch_gene_info_tsv            = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Report generation
    ------------------------------------------------------------------------------*/
    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    // Prepare MultiQC inputs
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    //
    // MODULE: MultiQC
    //
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    ch_multiqc_report  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html


    emit:
    // virus_enrichment_tsv        = ch_virus_enrichment_tsv
    // virus_classification_tsv    = ch_virus_summaries_tsv
    // virus_extension_tsv         = ch_virus_extension_summary_tsv // Don't want to update all other workflow snapshots
    // virus_quality_tsv           = ch_quality_summary_tsv
    // filtered_viruses_fna_gz     = ch_filtered_viruses_fna_gz
    // anicluster_reps_fna_gz      = ch_anicluster_reps_fasta_gz
    // alignment_results_tsv       = ch_alignment_results_tsv   // Inconsistent hash
    // host_predictions_tsv        = ch_host_predictions_tsv
    // marker_taxonomy_tsv         = ch_marker_taxonomy_tsv
    // bacphlip_lifestyle_tsv      = ch_bacphlip_lifestyle_tsv  // Inconsistent hash
    // pharokka_output_tsv         = ch_pharokka_output_tsv
    // instrain_gene_info          = ch_gene_info_tsv
    // propagate_results_tsv       = ch_provirus_activity_tsv
    multiqc_report              = ch_multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
