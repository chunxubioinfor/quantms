process PROTEOMICSLFQ_DETECT {
    tag "${meta.mzml_id}"
    label 'process_medium'
    label 'openms'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/bigbio/openms-tools-thirdparty-sif:2026.10.04' :
        'ghcr.io/bigbio/openms-tools-thirdparty:2026.10.04' }"

    input:
    tuple val(meta), path(mzml), path(id_file)
    path(expdes)
    path(fasta)

    output:
    path "feat_ckpt/*.featureParquet", type: 'dir', emit: checkpoint
    path "*.log", emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def decoys_present = params.quantify_decoys ? '-PeptideQuantification:quantify_decoys' : ''
    def feature_with_id_min_score =  "-feature_with_id_min_score ${params.feature_with_id_min_score}"
    def feature_without_id_min_score = params.targeted_only == false ? "-feature_without_id_min_score ${params.feature_without_id_min_score}" : ""

    // Per-run half of ProteomicsLFQ (OpenMS >= 3.6.0): detects features for this one run and writes
    // a checkpoint named after its design row; PROTEOMICSLFQ then combines all of them via -feat_dir.
    // The full design is passed on purpose -- the checkpoint records this run's design row from it.
    // Only detection-relevant options are passed here; they must match those PROTEOMICSLFQ passes,
    // otherwise the combine step refuses the checkpoint and names the differing setting.
    // Combine-only options (FDR, inference, protein quantification, alignment order) are left out so
    // that changing them reruns only the combine step on -resume.
    """
    ProteomicsLFQ \\
        -threads ${task.cpus} \\
        -in ${mzml} \\
        -ids ${id_file} \\
        -design ${expdes} \\
        -fasta ${fasta} \\
        -feat_dir feat_ckpt \\
        -detect_only \\
        -quantification_method ${params.quantification_method} \\
        -targeted_only ${params.targeted_only} \\
        ${feature_with_id_min_score} \\
        ${feature_without_id_min_score} \\
        -mass_recalibration ${params.mass_recalibration} \\
        -Seeding:algorithm ${params.lfq_seeding_algorithm} \\
        -Seeding:intThreshold ${params.lfq_intensity_threshold} \\
        ${decoys_present} \\
        $args \\
        2>&1 | tee ${meta.mzml_id}_proteomicslfq_detect.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ProteomicsLFQ: \$(ProteomicsLFQ 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
    END_VERSIONS
    """
}
