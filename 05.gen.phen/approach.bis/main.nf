nextflow.enable.dsl = 2

process CAASTOOLS_DISCOVERY {
    tag "${gene_id} | ${cfg_id}"

    stageInMode 'copy'
    publishDir path: { "${params.results_root}/${params.run_id}/caas/by-config/${cfg_id}" },
               mode: 'copy', overwrite: true

    input:
    tuple val(gene_id), path(alignment), val(cfg_id), path(config), path(caastools_dir)

    output:
    tuple val(cfg_id), path("${gene_id}.${cfg_id}.caas.tsv"), emit: results

    script:
    """
    "${params.python_command}" "${caastools_dir}/ct" discovery \
        -a "${alignment}" \
        -t "${config}" \
        -o "${gene_id}.${cfg_id}.caas.tsv" \
        --fmt "${params.caas_alignment_format}" \
        --patterns "${params.caas_patterns}" \
        --max_bg_gaps "${params.caas_max_bg_gaps}" \
        --max_fg_gaps "${params.caas_max_fg_gaps}" \
        --max_gaps "${params.caas_max_gaps}" \
        --max_gaps_per_position "${params.caas_max_gaps_per_position}" \
        --max_bg_miss "${params.caas_max_bg_miss}" \
        --max_fg_miss "${params.caas_max_fg_miss}" \
        --max_miss "${params.caas_max_miss}"

    # CAAStools creates no file when a gene has no CAAS. Emit a header-only
    # table so every scheduled gene/config pair has a traceable result.
    if [[ ! -s "${gene_id}.${cfg_id}.caas.tsv" ]]; then
        printf 'Gene\tTrait\tPosition\tSubstitution\tPvalue\tPattern\tFFGN\tFBGN\tGFG\tGBG\tMFG\tMBG\tFFG\tFBG\tMS\n' \
            > "${gene_id}.${cfg_id}.caas.tsv"
    fi
    """
}

process MERGE_CONFIG_RESULTS {
    tag "${cfg_id} | all genes"

    publishDir path: { "${params.results_root}/${params.run_id}/merged/by-config" },
               mode: 'copy', overwrite: true

    input:
    tuple val(cfg_id), path(caas_results)
    path merge_script

    output:
    tuple val(cfg_id), path("${cfg_id}.all-genes.caas.tsv"), emit: results

    script:
    """
    "${params.python_command}" "${merge_script}" \
        --input-dir . \
        --suffix '.caas.tsv' \
        --merged "${cfg_id}.all-genes.caas.tsv"
    """
}

process MERGE_AND_FILTER_RESULTS {
    tag 'all configs | merge and significance filter'

    publishDir path: { "${params.results_root}/${params.run_id}/merged" },
               mode: 'copy', overwrite: true

    input:
    path per_config_results
    path merge_script

    output:
    path 'caas.all-configs.all-genes.tsv', emit: merged
    path 'caas.all-configs.all-genes.significant.tsv', emit: significant

    script:
    """
    "${params.python_command}" "${merge_script}" \
        --input-dir . \
        --suffix '.all-genes.caas.tsv' \
        --merged 'caas.all-configs.all-genes.tsv' \
        --significant 'caas.all-configs.all-genes.significant.tsv' \
        --alpha "${params.significance_threshold}"
    """
}

workflow {
    if (!params.run_id) {
        error 'Missing --run_id. Use run_pipeline.sh or submit_pipeline_slurm.sh.'
    }
    if (!params.alignments) {
        error 'params.alignments is empty. Edit conf/cluster.config.'
    }
    if (!params.caas_config_glob) {
        error 'params.caas_config_glob is empty. Edit conf/cluster.config.'
    }

    caastools_dir = file(params.caastools_dir, checkIfExists: true)
    merge_script = file(params.merge_script, checkIfExists: true)

    alignments_ch = Channel
        .fromPath(params.alignments, checkIfExists: true)
        .filter { alignment -> !alignment.isDirectory() }
        .map { alignment ->
            tuple(alignment.name.replaceFirst(/\.[^.]+$/, ''), alignment)
        }

    caas_configs_ch = Channel
        .fromPath(params.caas_config_glob, checkIfExists: true)
        .filter { config -> !config.isDirectory() }
        .map { config ->
            tuple(config.name.replaceFirst(/\.caas\.cfg$/, ''), config)
        }

    caas_jobs_ch = alignments_ch
        .combine(caas_configs_ch)
        .map { gene_id, alignment, cfg_id, config ->
            tuple(gene_id, alignment, cfg_id, config, caastools_dir)
        }

    CAASTOOLS_DISCOVERY(caas_jobs_ch)

    // groupTuple without a fixed size emits only after the discovery channel
    // closes, so merging starts after every gene/config job has completed.
    results_by_config_ch = CAASTOOLS_DISCOVERY.out.results.groupTuple(by: 0)
    MERGE_CONFIG_RESULTS(results_by_config_ch, merge_script)

    per_config_results_ch = MERGE_CONFIG_RESULTS.out.results
        .map { cfg_id, result -> result }
        .collect()

    MERGE_AND_FILTER_RESULTS(per_config_results_ch, merge_script)
}
