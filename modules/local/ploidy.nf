// By João Pitta (jlpitta82@gmail.com) and Beatriz Toscano (beatriz.melo@fiocruz.br)
// At Fiocruz-PE
process PLOIDY_CHECK {
    tag { sample }
    label 'process_medium'
    publishDir { "${params.outdir}/${sample}/qc/ploidy" }, mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}.ploidy_call.json"), emit: call
    path "gs2_out", emit: genomescope_report
    path "${sample}_smudgeplot_report.json", emit: smudgeplot_report
    path "${sample}*.png", emit: plots

    script:
    // Delega pro ploidycheck externo (github.com/jlpitta/ploidycheck) em vez de
    // manter uma cópia embutida da lógica FastK→GenomeScope2→Smudgeplot→
    // ploidy_call.py — o fungiflow só usa o resultado (decisão de parâmetro do
    // montador), o método/critério de decisão (AB dominante + cobertura
    // batendo com kmercov, L=7) mora e é documentado no repo do ploidycheck.
    // O script gerencia o próprio ambiente conda/mamba internamente, por isso
    // não usa a diretiva `conda` do Nextflow aqui.
    def L = params.ploidy_smudge_l ?: 7
    def AB = params.ploidy_ab_threshold ?: 0.5
    def COV_TOL = params.ploidy_coverage_tolerance ?: 0.3
    """
    ${params.ploidycheck_bin} \
        --sample ${sample} \
        --r1 ${r1} \
        --r2 ${r2} \
        --outdir . \
        --threads ${task.cpus} \
        --smudge-l ${L} \
        --ab-fraction-threshold ${AB} \
        --coverage-tolerance ${COV_TOL} \
        --env ${params.ploidycheck_env ?: 'ploidycheck'}
    """
}
