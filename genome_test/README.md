# Dados de teste do fungiflow

Datasets pequenos para validar o pipeline de ponta a ponta sem precisar de dados próprios. Os dois primeiros são herdados do bacflow (genomas bacterianos, úteis apenas enquanto os módulos ainda não foram trocados pelos equivalentes fúngicos); os dois seguintes são datasets fúngicos do projeto — um haploide/homozigoto, outro heterozigótico sintético (par de haplótipos divergentes) para validar o branch de detecção de ploidia. Todos cabem em `--t 4` e rodam em poucos minutos numa máquina comum.

## `mycoplasma_genitalium_synthetic/` — rápido, sintético a partir de sequência real

- **Referência real:** *Mycoplasmoides genitalium* G37, [GCF_000027325.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000027325.1/) (580 kb — um dos menores genomas bacterianos conhecidos).
- **Reads:** simulados a partir dessa sequência real (não são reads de sequenciador de verdade) — ~40x long reads (erro aleatório ~5%, perfil simplificado, não replica o padrão sistemático real de erro ONT/PacBio) e ~30x short reads (erro ~0.5%).
- **Uso:** smoke test rápido — confirma que o pipeline roda de ponta a ponta e gera todos os relatórios de QC esperados. **Não é adequado para avaliar a qualidade real do polimento** — como os reads não têm o padrão de erro sistemático de um sequenciador de verdade, o Medaka já corrige praticamente tudo sozinho, então os relatórios pré/pós-polish tendem a sair idênticos (não é bug do pipeline, é limitação do dado sintético).

```bash
nextflow run fungiflow.nf --t 4 \
    --long_reads genome_test/mycoplasma_genitalium_synthetic/long_reads.fastq.gz \
    --short_reads_1 genome_test/mycoplasma_genitalium_synthetic/short_reads_1.fastq.gz \
    --short_reads_2 genome_test/mycoplasma_genitalium_synthetic/short_reads_2.fastq.gz \
    --genome_size 580000 \
    --sample_name mgenitalium_test \
    --reference genome_test/mycoplasma_genitalium_synthetic/reference.fasta
```

## `staphylococcus_aureus_real/` — mais lento, 100% dados reais de sequenciador

- **Long + short reads:** **mesma cepa real** (*Staphylococcus aureus* JH62PP1, amostra [SAMD00828832](https://www.ebi.ac.uk/ena/browser/view/SAMD00828832) no ENA/DDBJ) — long reads ONT ([DRR613158](https://www.ebi.ac.uk/ena/browser/view/DRR613158)) e short reads Illumina ([DRR613151](https://www.ebi.ac.uk/ena/browser/view/DRR613151)), ambos subamostrados com `seqkit sample` (seed 42) para ~25x de cobertura cada, a partir dos runs completos públicos.
- **Referência:** *Staphylococcus aureus* NCTC 8325, [GCF_000013425.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013425.1/) — **cepa diferente** da dos reads (de propósito: no uso real do pipeline raramente se tem uma referência da cepa exata, então isso testa o cenário realista de comparar contra uma referência próxima, não idêntica).
- **Uso:** validação de verdade do polimento com erro real de sequenciador. Resultado observado num run de teste (`--reference` informado, caminho Flye):

  | Métrica | Pré-polish | Pós-polish |
  |---|---|---|
  | Indels /100kbp (QUAST) | 119.96 | 55.69 |
  | CheckM2 Completeness | 90.7% | **100.0%** |
  | CheckM2 Contamination | 8.53% | **0.05%** |

  A queda de completude/contaminação pré-polish é o efeito esperado de indels de long-read causando frameshift (genes fragmentados); o polimento com short reads corrige isso quase totalmente. Mismatches (~1100/100kbp) refletem principalmente divergência genômica real entre as duas cepas, não erro de sequenciamento — por isso ficam altos mesmo pós-polish.

```bash
nextflow run fungiflow.nf --t 4 \
    --long_reads genome_test/staphylococcus_aureus_real/long_reads.fastq.gz \
    --short_reads_1 genome_test/staphylococcus_aureus_real/short_reads_1.fastq.gz \
    --short_reads_2 genome_test/staphylococcus_aureus_real/short_reads_2.fastq.gz \
    --genome_size 2.8m \
    --sample_name saureus_test \
    --reference genome_test/staphylococcus_aureus_real/reference.fasta
```

## `saccharomyces_cerevisiae_synthetic/` — primeiro dataset fúngico, sintético a partir de sequência real

- **Referência real:** *Saccharomyces cerevisiae* S288C, [GCF_000146045.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/) (montagem R64, ~12.16 Mb, 17 sequências — 16 cromossomos + mitocôndria). Genoma haploide "gold standard", pouco repetitivo — bom primeiro caso pra validar o pipeline sem a complexidade extra de um fungo filamentoso.
- **Reads:** simulados a partir dessa sequência real com [Badread](https://github.com/rrwick/Badread) (long reads, perfil `random`/`ideal` — erro simplificado, não replica o padrão sistemático real de ONT/PacBio) e `wgsim` (short reads, erro ~0.5%, seed 42).
  - Long reads: gerados a ~40x e depois subamostrados (`seqkit sample -p 0.175`, seed 42) para ~7x — o arquivo a 40x direto do Badread deu 469 MB, acima do limite de 100 MB por arquivo do GitHub.
  - Short reads: ~30x (150 bp pareado), sem subamostragem (78 MB por arquivo, dentro do limite).
- **Status dos módulos:** este dataset já existe, mas os módulos do pipeline (`unicycler.nf`, `checkm2.nf`, `bakta.nf`, `gtdbtk.nf`, `amrfinder.nf`) ainda são os originais do bacflow (bacterianos) — a Etapa C do plano (trocar por SPAdes/MaSuRCA, EukCC, Funannotate, ITSx+UNITE) ainda não foi feita. Rodar este dataset hoje só exercita o caminho de QC/montagem/polimento genéricos (Flye, Racon, Medaka, Polypolish, NextPolish, QUAST, BUSCO com lineage bacteriana ainda), não os módulos fúngico-específicos.
- **Uso:** smoke test rápido assim que os módulos da Etapa C forem trocados — mesmo papel que o `mycoplasma_genitalium_synthetic/` tem pro bacflow. **Não é adequado para avaliar qualidade real de polimento/anotação** pela mesma razão do dataset de Mycoplasma (erro sintético não replica o padrão real de sequenciador) e pela baixa cobertura de long reads (~7x).

```bash
nextflow run fungiflow.nf --t 4 \
    --long_reads genome_test/saccharomyces_cerevisiae_synthetic/long_reads.fastq.gz \
    --short_reads_1 genome_test/saccharomyces_cerevisiae_synthetic/short_reads_1.fastq.gz \
    --short_reads_2 genome_test/saccharomyces_cerevisiae_synthetic/short_reads_2.fastq.gz \
    --genome_size 12.16m \
    --sample_name scerevisiae_test \
    --reference genome_test/saccharomyces_cerevisiae_synthetic/reference.fasta
```

## `saccharomyces_cerevisiae_heterozygous/` — dataset heterozigótico/poliploide sintético

- **Objetivo:** validar o branch de detecção de ploidia do fungiflow (GenomeScope2/Smudgeplot pré-montagem, purge_dups/nQuire pós-montagem — ver plano de poliploidia na página Notion do projeto). O dataset `saccharomyces_cerevisiae_synthetic/` é essencialmente haploide/homozigoto (cepa de laboratório S288C) e não serve para isso.
- **Como foi construído:** a partir da mesma referência S288C/R64 ([GCF_000146045.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/)) usada no dataset haploide, tratada como haplótipo A (`hapA.fasta`). Um haplótipo B divergente (`hapB.fasta`) foi gerado com `scripts/mutate_haplotype.py` (script novo deste projeto, SNPs pontuais aleatórios, sem indels), introduzindo **1.5% de divergência** (182.321 SNPs, seed 42) — dentro da faixa típica de heterozigose reportada para leveduras heterozigóticas/híbridas na literatura.
- **Reads:** gerados independentemente para cada haplótipo e depois combinados num único par de arquivos, simulando o que se obteria sequenciando uma célula real com dois haplótipos distintos (amostragem aleatória de qualquer uma das cópias):
  - Short reads (`wgsim`, 150 bp pareado, erro ~0.5%): ~15x de hapA (seed 42) + ~15x de hapB (seed 43) = **~30x combinado** — mesma escala do dataset haploide.
  - Long reads (`badread`, perfil `random`/`ideal`): ~3.5x de hapA (seed 42) + ~3.5x de hapB (seed 43) = **~7x combinado** — mesma escala do dataset haploide; gerados diretamente na cobertura alvo, sem precisar de subamostragem pós-hoc.
- **Arquivos extra:** `hapA.fasta`/`hapB.fasta` (os dois haplótipos "verdadeiros", mantidos para permitir conferir se a chamada de ploidia bate com a divergência real introduzida). `reference.fasta` é uma cópia de `hapA.fasta`, mantida só por convenção de compatibilidade com o pipeline (`--reference`); não representa uma referência "correta" única para um genoma heterozigótico de verdade.
- **Status:** dataset pronto (07/08/2026); ainda não validado por nenhuma ferramenta de ploidia, porque essas ferramentas (GenomeScope2, Smudgeplot, purge_dups, nQuire) ainda não foram implementadas no pipeline — esse dataset é o desbloqueador (passo 1 do plano de poliploidia), não uma validação em si.
- **Resultado esperado quando o branch de ploidia existir:** GenomeScope2/Smudgeplot devem reportar heterozigose ≈1.5% e um padrão de ploidia `AB` (diploide), não `AA` (haploide) como no dataset `_synthetic`.

```bash
nextflow run fungiflow.nf --t 4 \
    --long_reads genome_test/saccharomyces_cerevisiae_heterozygous/long_reads.fastq.gz \
    --short_reads_1 genome_test/saccharomyces_cerevisiae_heterozygous/short_reads_1.fastq.gz \
    --short_reads_2 genome_test/saccharomyces_cerevisiae_heterozygous/short_reads_2.fastq.gz \
    --genome_size 12.16m \
    --sample_name scerevisiae_het_test \
    --reference genome_test/saccharomyces_cerevisiae_heterozygous/reference.fasta
```

## Sem `--reference`

Todos os comandos acima também podem ser rodados sem `--reference` (removendo a flag) para exercitar o caminho BUSCO em vez da comparação QUAST baseada em referência.
