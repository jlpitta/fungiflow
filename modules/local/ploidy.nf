// By João Pitta (jlpitta82@gmail.com) and Beatriz Toscano (beatriz.melo@fiocruz.br)
// At Fiocruz-PE
process FASTK_HIST {
    tag { sample }
    label 'process_medium'
    conda "${System.getenv('HOME')}/miniforge3/envs/fungiflow-ploidy"
    publishDir { "${params.outdir}/${sample}/qc/ploidy" }, mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_fastk.ktab"), path("${sample}.histo"), emit: ktab_histo

    script:
    """
    FastK -k21 -T${task.cpus} -M8 -t1 -N${sample}_fastk ${r1} ${r2}
    Histex -A ${sample}_fastk > ${sample}.histo
    """
}

process GENOMESCOPE2 {
    tag { sample }
    label 'process_low'
    conda "${System.getenv('HOME')}/miniforge3/envs/fungiflow-ploidy"
    publishDir { "${params.outdir}/${sample}/qc/ploidy" }, mode: 'copy'

    input:
    tuple val(sample), path(ktab), path(histo)

    output:
    tuple val(sample), path("gs2_out"), emit: report

    script:
    """
    genomescope2 -i ${histo} -o gs2_out -k 21 -p 2 -n ${sample}
    """
}

process SMUDGEPLOT {
    tag { sample }
    label 'process_medium'
    conda "${System.getenv('HOME')}/miniforge3/envs/fungiflow-ploidy"
    publishDir { "${params.outdir}/${sample}/qc/ploidy" }, mode: 'copy'

    input:
    tuple val(sample), path(ktab), path(histo)

    output:
    tuple val(sample), path("${sample}_smudgeplot_report.json"), emit: report
    path "${sample}*.png", emit: plots

    script:
    // L=7: cutoff empírico validado contra o kmercov do GenomeScope2 nos
    // datasets de teste (haploide e heterozigótico, ~11-12x de cobertura de
    // k-mer) — ver seção 7 do documento "fungiflow" no Notion. O cutoff que o
    // próprio `smudgeplot cutoff` sugere para esses datasets (~10) infla a
    // cobertura 1n inferida em quase 2x mantendo uma smudge AB enganosamente
    // "limpa" — não é confiável nessa faixa de cobertura, por isso não é usado
    // como default aqui apesar de ser a recomendação da ferramenta.
    def L = params.ploidy_smudge_l ?: 7
    """
    smudgeplot hetmers -L ${L} -t ${task.cpus} -o ${sample}_kmerpairs ${ktab}
    smudgeplot all -o ${sample} -t ${sample} --json_report ${sample}_kmerpairs.smu
    """
}

process PLOIDY_CALL {
    tag { sample }
    label 'process_low'
    conda "${System.getenv('HOME')}/miniforge3/envs/fungiflow-ploidy"
    publishDir { "${params.outdir}/${sample}/qc/ploidy" }, mode: 'copy'

    input:
    tuple val(sample), path(gs2_report), path(smudge_report)

    output:
    tuple val(sample), path("${sample}.ploidy_call.json"), emit: call

    script:
    // heterozigose só é reportada quando os dois sinais concordam: a smudge AB
    // é maioria da massa (>= --ab-fraction-threshold) E a cobertura 1n que o
    // Smudgeplot infere bate com o kmercov do GenomeScope2 (dentro de
    // --coverage-tolerance) — ver docstring de ploidy_call.py.
    """
    ploidy_call.py \
        --sample ${sample} \
        --genomescope-model ${gs2_report}/${sample}_model.txt \
        --smudgeplot-json ${smudge_report} \
        --ab-fraction-threshold ${params.ploidy_ab_threshold ?: 0.5} \
        --coverage-tolerance ${params.ploidy_coverage_tolerance ?: 0.3} \
        --out ${sample}.ploidy_call.json
    """
}
