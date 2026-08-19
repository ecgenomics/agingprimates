nextflow.enable.dsl = 2

process PREPARE_POOLED_HYPOTHESES {
    tag "${params.pooled_comparisons}-of-max"

    stageInMode 'copy'
    publishDir path: { "${params.results_root}/${params.run_id}/metadata" },
               mode: 'copy', overwrite: true

    input:
    path pool_config
    path caastools_dir
    path preparation_script

    output:
    path "longevity.pooled.hypotheses.tsv", emit: hypotheses
    path "longevity.pooled.hypotheses.metadata.tsv", emit: metadata

    script:
    """
    "${params.python_command}" "${preparation_script}" \
        --caastools-dir "${caastools_dir}" \
        --pool-file "${pool_config}" \
        --output "longevity.pooled.hypotheses.tsv" \
        --metadata-output "longevity.pooled.hypotheses.metadata.tsv" \
        --fg-size "${params.pooled_fg_size}" \
        --bg-size "${params.pooled_bg_size}" \
        --comparisons "${params.pooled_comparisons}" \
        --seed "${params.pooled_seed}"
    """
}

process CAASTOOLS_POOLED_DISCOVERY {
    tag "${gene_id}"

    stageInMode 'copy'
    publishDir path: { "${params.results_root}/${params.run_id}/caas-pooled" },
               pattern: "*.pooled.caas.tsv", mode: 'copy', overwrite: true
    publishDir path: { "${params.results_root}/${params.run_id}/caas-pooled-events" },
               pattern: "*.pooled.caas.events.tsv", mode: 'copy', overwrite: true

    input:
    tuple val(gene_id), path(alignment)
    path pool_config
    path pooled_hypotheses
    path caastools_dir

    output:
    tuple val(gene_id), path("${gene_id}.pooled.caas.tsv"), emit: legacy_results
    tuple val(gene_id), path("${gene_id}.pooled.caas.events.tsv"), emit: event_results

    script:
    """
    "${params.python_command}" "${caastools_dir}/ct" pooled-discovery \
        -a "${alignment}" \
        -t "${pool_config}" \
        -s "${pooled_hypotheses}" \
        -o "${gene_id}.pooled.caas.tsv" \
        --event-output "${gene_id}.pooled.caas.events.tsv" \
        --hypotheses-output none \
        --fmt "${params.caas_alignment_format}" \
        --filter_significant "${params.caas_filter_significant}" \
        --patterns "${params.caas_patterns}" \
        --max_bg_gaps "${params.caas_max_bg_gaps}" \
        --max_fg_gaps "${params.caas_max_fg_gaps}" \
        --max_gaps "${params.caas_max_gaps}" \
        --max_gaps_per_position "${params.caas_max_gaps_per_position}" \
        --max_bg_miss "${params.caas_max_bg_miss}" \
        --max_fg_miss "${params.caas_max_fg_miss}" \
        --max_miss "${params.caas_max_miss}" \
        --min_fg_observed "${params.caas_min_fg_observed}" \
        --min_bg_observed "${params.caas_min_bg_observed}"
    """
}

workflow {
    if (!params.run_id) {
        error "Missing --run_id. Use run_pipeline.sh or submit_pipeline_slurm.sh."
    }
    if (!params.alignments) {
        error "params.alignments is empty. Edit conf/cluster.config."
    }
    if (!params.pooled_trait_config) {
        error "params.pooled_trait_config is empty. Edit conf/cluster.config."
    }

    pool_config = file(params.pooled_trait_config, checkIfExists: true)
    caastools_dir = file(params.caastools_dir, checkIfExists: true)
    preparation_script = file("${projectDir}/scripts/prepare_pooled_hypotheses.py", checkIfExists: true)

    pooled = PREPARE_POOLED_HYPOTHESES(pool_config, caastools_dir, preparation_script)
    pooled_hypotheses = pooled.hypotheses

    pooled_jobs_ch = Channel
        .fromPath(params.alignments, checkIfExists: true)
        .filter { alignment -> !alignment.isDirectory() }
        .map { alignment ->
            tuple(
                alignment.name.replaceFirst(/\.[^.]+$/, ''),
                alignment
            )
        }

    CAASTOOLS_POOLED_DISCOVERY(
        pooled_jobs_ch,
        pool_config,
        pooled_hypotheses,
        caastools_dir
    )
}
