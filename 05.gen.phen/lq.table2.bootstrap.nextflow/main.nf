nextflow.enable.dsl = 2

process CAASTOOLS_BOOTSTRAP {
    tag "${gene_id}"

    stageInMode 'copy'
    publishDir path: { "${params.results_root}/${params.run_id}/caas-bootstrap" },
               mode: 'copy', overwrite: true

    input:
    tuple val(gene_id), path(alignment), path(trait_config), path(resampled_traits), path(caastools_dir)

    output:
    tuple val(gene_id), path("${gene_id}.bootstrap.caas.tsv"), emit: results

    script:
    """
    "${params.python_command}" "${caastools_dir}/ct" bootstrap \
        -a "${alignment}" \
        -t "${trait_config}" \
        -s "${resampled_traits}" \
        -o "${gene_id}.bootstrap.caas.tsv" \
        --fmt "${params.caas_alignment_format}" \
        --filter_significant "${params.caas_filter_significant}" \
        --patterns "${params.caas_patterns}" \
        --max_bg_gaps "${params.caas_max_bg_gaps}" \
        --max_fg_gaps "${params.caas_max_fg_gaps}" \
        --max_gaps "${params.caas_max_gaps}" \
        --max_gaps_per_position "${params.caas_max_gaps_per_position}" \
        --max_bg_miss "${params.caas_max_bg_miss}" \
        --max_fg_miss "${params.caas_max_fg_miss}" \
        --max_miss "${params.caas_max_miss}"
    """
}

workflow {
    if (!params.run_id) {
        error "Missing --run_id. Use run_pipeline.sh or submit_pipeline_slurm.sh."
    }
    if (!params.alignments) {
        error "params.alignments is empty. Edit conf/cluster.config."
    }
    if (!params.bootstrap_trait_config) {
        error "params.bootstrap_trait_config is empty. Edit conf/cluster.config."
    }
    if (!params.bootstrap_resampled_traits) {
        error "params.bootstrap_resampled_traits is empty. Edit conf/cluster.config."
    }

    trait_config = file(params.bootstrap_trait_config, checkIfExists: true)
    resampled_traits = file(params.bootstrap_resampled_traits, checkIfExists: true)
    caastools_dir = file(params.caastools_dir, checkIfExists: true)

    bootstrap_jobs_ch = Channel
        .fromPath(params.alignments, checkIfExists: true)
        .filter { alignment -> !alignment.isDirectory() }
        .map { alignment ->
            tuple(
                alignment.name.replaceFirst(/\.[^.]+$/, ''),
                alignment,
                trait_config,
                resampled_traits,
                caastools_dir
            )
        }

    CAASTOOLS_BOOTSTRAP(bootstrap_jobs_ch)
}
