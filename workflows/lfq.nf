/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES: Local to the pipeline
//
include { PROTEOMICSLFQ } from '../modules/local/openms/proteomicslfq/main'
include { PROTEOMICSLFQ_DETECT } from '../modules/local/openms/proteomicslfq_detect/main'
include { QPX_OPENMSCONSENSUS } from '../modules/bigbio/qpx/openmsconsensus/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { ID } from '../subworkflows/local/id/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/


workflow LFQ {
    take:
    ch_file_preparation_results
    ch_expdesign
    ch_database_wdecoy

    main:

    ch_software_versions = channel.empty()

    //
    // SUBWORKFLOWS: ID
    //
    ID(ch_file_preparation_results, ch_database_wdecoy, ch_expdesign)
    ch_software_versions = ch_software_versions.mix(ID.out.versions)

    //
    // SUBWORKFLOW: PROTEOMICSLFQ
    //
    ch_mzml_ids = ch_file_preparation_results.join(ID.out.id_results)

    if (params.plfq_parallel_detection && params.quantification_method == 'spectral_counting') {
        error("--plfq_parallel_detection requires --quantification_method feature_intensity: spectral counting detects no features.")
    }

    if (params.plfq_parallel_detection) {
        //
        // MODULE: PROTEOMICSLFQ_DETECT -- one feature-detection task per run, writing a
        // checkpoint each; PROTEOMICSLFQ then combines all of them (no mzML/IDs needed).
        // The combine can only start once every run of the design has a checkpoint, hence collect().
        //
        PROTEOMICSLFQ_DETECT(ch_mzml_ids, ch_expdesign.first(), ch_database_wdecoy.first())
        ch_software_versions = ch_software_versions.mix(PROTEOMICSLFQ_DETECT.out.versions.first())

        PROTEOMICSLFQ([],
                    [],
                    ch_expdesign,
                    ch_database_wdecoy,
                    PROTEOMICSLFQ_DETECT.out.checkpoint.collect()
                )
    } else {
        ch_mzml_ids
            .multiMap { it ->
                mzmls: pmultiqc_mzmls: it[1]
                ids: it[2]
            }
            .set{ ch_plfq }
        PROTEOMICSLFQ(ch_plfq.mzmls.collect(),
                    ch_plfq.ids.collect(),
                    ch_expdesign,
                    ch_database_wdecoy,
                    []
                )
    }
    ch_software_versions = ch_software_versions.mix(PROTEOMICSLFQ.out.versions)

    //
    // MODULE: QPX_OPENMSCONSENSUS  -- convert the OpenMS consensusXML + SDRF into the
    // final clean QPX dataset + MuData (the published quantification artifact).
    //
    QPX_OPENMSCONSENSUS(
        PROTEOMICSLFQ.out.out_consensusXML,
        file(params.input),
        params.accession ?: '',
        // The database the search used. qpx fills null pg.sequence_coverage,
        // pg.molecular_weight and feature.pg_positions from it; decoy entries are
        // skipped, so the with-decoy database is safe to pass. .first() keeps this
        // a value channel, since the same channel already feeds ID and PROTEOMICSLFQ.
        ch_database_wdecoy.first(),
    )
    ch_software_versions = ch_software_versions.mix(QPX_OPENMSCONSENSUS.out.versions)

    ID.out.psmrescoring_results
        .map { it -> it[1] }
        .set { ch_pmultiqc_ids }

    ID.out.ch_consensus_results
        .map { it -> it[1] }
        .set { ch_pmultiqc_consensus }

    emit:
    ch_pmultiqc_ids         = ch_pmultiqc_ids
    ch_pmultiqc_consensus   = ch_pmultiqc_consensus
    final_result            = QPX_OPENMSCONSENSUS.out.qpx_dataset
    versions                = ch_software_versions
    msstats_in              = PROTEOMICSLFQ.out.out_msstats
}

/*
========================================================================================
    THE END
========================================================================================
*/
