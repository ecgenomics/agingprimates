nextflow.enable.dsl = 2

process CAASTOOLS_DISCOVERY {
    tag "${gene_id} | ${cfg_id}"

    stageInMode 'copy'
    publishDir path: { "${params.results_root}/${params.run_id}/caas/${cfg_id}" },
               mode: 'copy', overwrite: true

    input:
    tuple val(gene_id), path(alignment), val(cfg_id), path(config), path(caastools_dir)

    output:
    tuple val(gene_id), val(cfg_id), path("${gene_id}.${cfg_id}.caas"), emit: results

    script:
    """
    "${params.python_command}" "${caastools_dir}/ct" discovery \
        -a "${alignment}" \
        -t "${config}" \
        -o "${gene_id}.${cfg_id}.caas" \
        --fmt "${params.caas_alignment_format}"

    if [[ ! -f "${gene_id}.${cfg_id}.caas" ]]; then
        touch "${gene_id}.${cfg_id}.caas"
    fi
    """
}

process RERCONVERGE {
    tag "${cfg_id}"

    stageInMode 'copy'
    publishDir path: { "${params.results_root}/${params.run_id}/rerconverge/${cfg_id}" },
               mode: 'copy', overwrite: true

    input:
    tuple val(cfg_id), path(config), path(tree_manifest), path(master_tree), path(rer_script)

    output:
    tuple val(cfg_id), path("rerconverge.${cfg_id}"), emit: results

    script:
    """
    "${params.rscript_command}" "${rer_script}" \
        --trees "${tree_manifest}" \
        --master-tree "${master_tree}" \
        --phenotype "${config}" \
        --outdir "rerconverge.${cfg_id}" \
        --clade "${params.rer_clade}" \
        --transition "${params.rer_transition}" \
        --weighted "${params.rer_weighted}" \
        --transform "${params.rer_transform}" \
        --impute "${params.rer_impute}" \
        --min-trees "${params.rer_min_trees}" \
        --max-trees "${params.rer_max_trees}" \
        --min-species "${params.rer_min_species}" \
        --min-valid "${params.rer_min_valid}" \
        --min-foreground "${params.rer_min_foreground}" \
        --bootstrap "${params.rer_bootstrap}" \
        --bootn "${params.rer_bootn}"
    """
}

workflow {
    if (!params.run_id) {
        error "Missing --run_id. Use run_pipeline.sh or submit_pipeline_slurm.sh."
    }
    if (!params.run_caas && !params.run_rerconverge) {
        error "Both run_caas and run_rerconverge are false; there is nothing to run."
    }

    if (params.run_caas) {
        if (!params.alignments) {
            error "run_caas=true but params.alignments is empty. Edit conf/cluster.config."
        }
        if (!params.caas_config_glob) {
            error "run_caas=true but params.caas_config_glob is empty."
        }

        caastools_dir = file(params.caastools_dir, checkIfExists: true)

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
    }

    if (params.run_rerconverge) {
        if (!params.rerconverge_config_glob) {
            error "run_rerconverge=true but params.rerconverge_config_glob is empty."
        }
        if (!params.rerconverge_trees) {
            error "run_rerconverge=true but params.rerconverge_trees is empty. Edit conf/cluster.config."
        }
        if (!params.rerconverge_master_tree) {
            error "run_rerconverge=true but params.rerconverge_master_tree is empty. Edit conf/cluster.config."
        }

        tree_manifest = file(params.rerconverge_trees, checkIfExists: true)
        master_tree = file(params.rerconverge_master_tree, checkIfExists: true)
        rer_script = file(params.rerconverge_script, checkIfExists: true)

        rer_jobs_ch = Channel
            .fromPath(params.rerconverge_config_glob, checkIfExists: true)
            .filter { config -> !config.isDirectory() }
            .map { config ->
                tuple(
                    config.name.replaceFirst(/\.rerc\.cfg$/, ''),
                    config,
                    tree_manifest,
                    master_tree,
                    rer_script
                )
            }

        RERCONVERGE(rer_jobs_ch)
    }
}
