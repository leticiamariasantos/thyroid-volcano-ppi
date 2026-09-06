# ═══════════════════════════════════════════════════════════════════════════════
# 04_panel.R — Definição a priori do painel de 30 vias (10 originais + 20 novas)
#
# A seleção das 20 novas vias é FEITA ANTES de qualquer interpretação dos
# resultados desta execução, com base em critérios biológicos objetivos:
#   (1) relevância p/ carcinoma de tireoide; (2) relevância p/ biologia tumoral;
#   (3) complementaridade às 10 vias originais; (4) mecanismos celulares
#   distintos; (5) qualidade/estabilidade da definição; (6) disponibilidade
#   confiável (KEGG/Reactome); (7) baixa redundância excessiva.
# Data da definição do painel: 2026-09-06 (mesma da nova execução).
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(KEGGREST)
  library(data.table)
})

log_msg("══ Definição do painel de 30 vias (a priori) ══")

# ───────────────────────────────────────────────────────────────────────────────
# Painel ORIGINAL de 10 vias (INALTERADO — pré-especificado no projeto)
# ───────────────────────────────────────────────────────────────────────────────
original10 <- data.table(
  panel_group = "ORIGINAL_10",
  database_id = c("hsa05216","hsa04919","hsa04010","hsa04151","hsa04150",
                  "hsa04115","hsa04210","hsa04110","hsa04310","hsa04064"),
  pathway_name = c("Thyroid cancer","Thyroid hormone signaling","MAPK signaling pathway",
                   "PI3K-Akt signaling pathway","mTOR signaling pathway","p53 signaling pathway",
                   "Apoptosis","Cell cycle","Wnt signaling pathway","NF-kappa B signaling pathway"),
  biological_domain = c("Oncogenesis (thyroid)","Hormone signaling","Growth factor signaling",
                        "Growth factor signaling","Nutrient signaling","Tumor suppression",
                        "Cell death","Cell proliferation","Developmental signaling","Inflammation/immune"),
  reason_for_inclusion = "Pré-especificado no desenho original do projeto (painel a priori)",
  redundancy_note = "—"
)

# ───────────────────────────────────────────────────────────────────────────────
# 20 vias ADICIONAIS — selecionadas a priori em 2026-09-06 (antes da análise)
# ───────────────────────────────────────────────────────────────────────────────
additional20 <- data.table(
  panel_group = "ADDITIONAL_20",
  database_id = c(
    "hsa04510",          # Focal adhesion
    "hsa04512",          # ECM-receptor interaction
    "R-HSA-1474244",     # Extracellular matrix organization
    "hsa04350",          # TGF-beta signaling
    "hsa04390",          # Hippo signaling
    "hsa04630",          # JAK-STAT signaling
    "hsa04668",          # TNF signaling
    "hsa04657",          # IL-17 signaling
    "hsa04612",          # Antigen processing and presentation
    "R-HSA-913531",      # Interferon Signaling
    "hsa04620",          # Toll-like receptor signaling
    "hsa04216",          # Ferroptosis
    "hsa04218",          # Cellular senescence
    "hsa03050",          # Proteasome
    "hsa00190",          # Oxidative phosphorylation
    "hsa03440",          # Homologous recombination
    "hsa03030",          # DNA replication
    "R-HSA-3299685",     # Detoxification of Reactive Oxygen Species
    "hsa04066",          # HIF-1 signaling
    "hsa04140"           # Autophagy - animal
  ),
  pathway_name = c(
    "Focal adhesion","ECM-receptor interaction","Extracellular matrix organization",
    "TGF-beta signaling pathway","Hippo signaling pathway","JAK-STAT signaling pathway",
    "TNF signaling pathway","IL-17 signaling pathway","Antigen processing and presentation",
    "Interferon Signaling","Toll-like receptor signaling pathway","Ferroptosis",
    "Cellular senescence","Proteasome","Oxidative phosphorylation",
    "Homologous recombination","DNA replication",
    "Detoxification of Reactive Oxygen Species","HIF-1 signaling pathway","Autophagy - animal"
  ),
  biological_domain = c(
    "Cell-matrix / adhesion","ECM / adhesion","ECM / collagen remodeling",
    "EMT / TGF-beta","Tissue growth control","Cytokine signaling",
    "Inflammation / TNF","Inflammation / IL-17","Immune / antigen presentation",
    "Immune / interferon","Innate immunity","Regulated cell death (iron)",
    "Cell cycle arrest / aging","Protein degradation","Energy metabolism / OXPHOS",
    "DNA repair","DNA replication","Oxidative stress / ROS",
    "Hypoxia / angiogenesis","Autophagy / stress response"
  ),
  reason_for_inclusion = c(
    "Medeia adesão célula-matriz e sinalização por integrinas; eixo central de invasão/migração tumoral",
    "Receptores de matriz (integrinas/colágeno) governam interação tumor-microambiente",
    "Biossíntese/remodelamento de colágeno e matriz; complementa as vias de sinalização de adesão com o componente estrutural",
    "Indutor canônico de EMT e imunossupressão no microambiente tumoral",
    "Regula tamanho de órgão e supressão tumoral (YAP/TAZ); mecanismo distinto das vias de crescimento originais",
    "Transdução de citocinas com papel em inflamação e imunidade antitumoral",
    "Inflamação crônica e sobrevivência/apoptose mediada por TNF no microambiente",
    "Eixo pró-inflamatório distinto do TNF, relevante para inflamação tecidual",
    "Mecanismo de apresentação de antígenos; conecta imunidade adaptativa ao reconhecimento tumoral",
    "Programa antiviral/imunológico de interferon, modulador do microambiente imune",
    "Sensoriamento imune inato (PAMPs); relevante para inflamação associada a tumor",
    "Morte celular dependente de ferro; vulnerabilidade metabólica distinta de apoptosis",
    "Parada proliferativa irreversível associada a senescência tumoral e resposta a terapia",
    "Degradação proteica ubiquitina-proteassoma; regula turnover de oncoproteínas",
    "Metabolismo energético mitocondrial; reprogramação metabólica tumoral (efeito Warburg inverso)",
    "Reparo de quebras de fita dupla; instabilidade genômica e resposta a dano de DNA",
    "Replicação do DNA; proliferação e estresse replicativo",
    "Eliminação de espécies reativas de oxigênio; homeostase redox",
    "Resposta a hipóxia e angiogênese; adaptação tumoral a microambiente hipóxico",
    "Reciclagem celular e sobrevivência sob estresse metabólico; mecanismo distinto de morte celular"
  ),
  redundancy_note = c(
    "Relacionada a ECM-receptor e ECM organization, mas representa o braço de sinalização",
    "Relacionada a Focal adhesion, mas representa o braço de receptores estruturais",
    "Relacionada a Focal adhesion/ECM-receptor, mas representa remodelamento estrutural de colágeno",
    "Distinta das vias de crescimento; converge com EMT",
    "Complementar a Wnt (original); eixo YAP/TAZ distinto",
    "Distinta de MAPK/PI3K; via própria de citocinas",
    "Distinta de NF-kB original embora relacionada; foco em TNF",
    "Distinta de TNF e NF-kB; eixo Th17",
    "Distinta das vias de sinalização; foco em MHC/processamento",
    "Distinta de TNF/TLR; programa de interferon tipo I",
    "Distinta de interferon; sensoriamento inato de PAMPs",
    "Morte celular não-apoptótica; distinta de Apoptosis original",
    "Distinta de Cell cycle/Apoptosis; parada irreversível",
    "Distinta de Apoptosis; degradação proteica",
    "Metabolismo mitocondrial; distinto das vias de sinalização",
    "Distinta de Cell cycle; reparo de DNA",
    "Distinta de Cell cycle; maquinaria replicativa",
    "Homeostase redox; distinta de OXPHOS",
    "Adaptação a hipóxia; distinta de OXPHOS e angiogênese",
    "Distinta de Apoptosis/Ferroptosis; reciclagem celular"
  )
)

# ── 1. Busca dos gene sets KEGG via KEGGREST ──────────────────────────────────
kegg_ids <- c(original10$database_id, additional20$database_id[grepl("^hsa", additional20$database_id)])

fetch_kegg_genes <- function(id) {
  res <- tryCatch(keggGet(id), error = function(e) NULL)
  if (is.null(res) || length(res) == 0) return(character(0))
  g <- res[[1]]$GENE
  if (is.null(g)) return(character(0))
  # GENE: posições ímpares = entrez; pares = "SYMBOL; desc"
  idx <- seq(2, length(g), by = 2)
  syms <- vapply(idx, function(i) {
    s <- strsplit(g[i], ";")[[1]][1]
    trimws(s)
  }, character(1))
  unique(syms[nzchar(syms)])
}

kegg_sets <- lapply(kegg_ids, fetch_kegg_genes)
names(kegg_sets) <- kegg_ids

# ── 2. Gene sets Reactome via GMT baixado ──────────────────────────────────────
reactome_gmt <- file.path(DIR_EXTERNAL, "reactome", "ReactomePathways.gmt")
gmt_lines <- readLines(reactome_gmt)
parse_gmt <- function(id) {
  hit <- grep(paste0("\t", id, "\t"), gmt_lines, fixed = TRUE)
  if (length(hit) == 0) return(character(0))
  parts <- strsplit(gmt_lines[hit[1]], "\t")[[1]]
  genes <- parts[-c(1, 2)]
  unique(genes[nzchar(genes)])
}
reactome_ids <- additional20$database_id[grepl("^R-HSA", additional20$database_id)]
reactome_sets <- lapply(reactome_ids, parse_gmt)
names(reactome_sets) <- reactome_ids

# ── 3. Montagem do painel final (30 vias) ─────────────────────────────────────
panel <- rbind(original10, additional20)
panel$source <- fifelse(grepl("^R-HSA", panel$database_id), "Reactome", "KEGG")

genesets <- c(kegg_sets, reactome_sets)
panel$geneset_size <- vapply(panel$database_id, function(id) length(genesets[[id]]), integer(1))
panel$geneset <- vapply(panel$database_id, function(id) paste(sort(genesets[[id]]), collapse = ","), character(1))

# validação
stopifnot(nrow(panel) == 30)
stopifnot(sum(panel$panel_group == "ORIGINAL_10") == 10)
stopifnot(sum(panel$panel_group == "ADDITIONAL_20") == 20)
stopifnot(all(panel$geneset_size > 0))
stopifnot(!anyDuplicated(panel$database_id))

fwrite_tsv(panel[, c("panel_group","source","database_id","pathway_name",
                     "biological_domain","reason_for_inclusion","redundancy_note","geneset_size")],
           file.path(DIR_PATH, "PANEL_30_PATHWAYS.tsv"))

# versão com a lista de genes completa (para auditoria/GSEA)
fwrite_tsv(panel, file.path(DIR_PATH, "PANEL_30_PATHWAYS_with_genes.tsv"))

# salva gene sets como lista nomeada (para fgsea) em formato leve
saveRDS(genesets, file.path(DIR_DATAOUT, "panel30_genesets.rds"))

log_msg("Painel final:", nrow(panel), "vias (10 originais + 20 adicionais)")
log_msg("  KEGG:", sum(panel$source == "KEGG"), "| Reactome:", sum(panel$source == "Reactome"))
log_msg("  Tamanho dos gene sets (min/med/max):",
        min(panel$geneset_size), median(panel$geneset_size), max(panel$geneset_size))
cat("\n")
print(panel[, .(database_id, pathway_name, source, geneset_size)])

log_msg("══ Painel definido e gene sets obtidos ══")
