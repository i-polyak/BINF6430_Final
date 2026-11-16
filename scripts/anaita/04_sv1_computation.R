# 04_sv1_computation.R
# BINF 6430 reproducibility project — recompute SV1 for the 33-sample GEO run
#
# Corrected to match the authors' only available SVA template (their human
# GTEx script, DifferentialGeneExpressionAnalysis_human.Rmd) as closely as
# possible, since their mouse SVA script was never shared:
#   - residualize expression against the condition model before EstDimRMT
#     (NOT just row-centering)
#   - explicit smartsva.cpp(B=1000, alpha=1), matching their call
#   - filter genes on rowSums(counts) > 1, matching their filter exactly
#   - seed set to 6430, per the group's Friday-meeting agreement
#
# NOTE: the authors' own script computes DESeq2 size factors but then takes
# log2() of the RAW (unnormalized) counts, not the size-factor-scaled matrix —
# the normalized matrix appears to be built and never used. This script
# replicates that literal behavior for a true reproduction. See the
# ALTERNATIVE block below if you want the presumably-intended version instead.
#
# Run as an sbatch job — see 04_submit_sv1.sbatch

suppressPackageStartupMessages({
  library(DESeq2)
  library(SmartSVA)
})

set.seed(6430)

## --- load tximport result from stage 1 ---
txi <- readRDS("~/binf6430_repro_env/txi.rds")
counts_raw <- round(txi$counts)  # DESeq2/sva-style tools expect integer-like counts

## --- Stage 1: build the 33-sample condition table ---
# Note: unlike Project 07's Reagan-data workbook, where the pairing/batch
# covariate (replicate: donor or cell-line identity) is KNOWN and simply
# added to the design formula directly, here the equivalent unwanted
# variation is UNKNOWN/latent -- that's why this project needs SV1 estimated
# via SmartSVA at all. If the Bloom authors had recorded and shared what
# SV1 actually represents (batch, dissection, RNA quality -- see the
# within-region SV1 pattern noted in Stage 6 below), this could have been a
# directly-modeled covariate like Project 07's `replicate`, with no SVA step
# needed.
coldata_33 <- data.frame(
  SRR = c("SRR32154025","SRR32154024","SRR32154023","SRR32154022","SRR32154021",
          "SRR32154020","SRR32154019","SRR32154018","SRR32154017","SRR32154016",
          "SRR32154015","SRR32154014","SRR32154013","SRR32154012","SRR32154011",
          "SRR32154010","SRR32154008","SRR32154007","SRR32154006","SRR32154005",
          "SRR32154004","SRR32154003","SRR32154002","SRR32154001","SRR32153999",
          "SRR32153998","SRR32153996","SRR32153995","SRR32153994","SRR32153993",
          "SRR32154009","SRR32154000","SRR32153997"),
  condition = c("PH_intact","PH_intact","PH_14dOVX","PH_14dOVX","PH_14dOVX",
                "PH_4mOVX","PH_4mOVX","PH_4mOVX","POA_intact","POA_intact",
                "POA_14dOVX","POA_14dOVX","POA_14dOVX","POA_4mOVX","POA_4mOVX",
                "POA_4mOVX","POA_intact","POA_intact","POA_14dOVX","POA_14dOVX",
                "POA_14dOVX","POA_4mOVX","POA_4mOVX","POA_4mOVX","PH_intact",
                "PH_14dOVX","PH_14dOVX","PH_4mOVX","PH_4mOVX","PH_4mOVX",
                "POA_intact","PH_intact","PH_14dOVX"),
  stringsAsFactors = FALSE
)
rownames(coldata_33) <- coldata_33$SRR
coldata_33 <- coldata_33[colnames(counts_raw), ]
stopifnot(all(rownames(coldata_33) == colnames(counts_raw)))
coldata_33$condition <- factor(coldata_33$condition)

cat("Condition table built. Counts per group:\n")
print(table(coldata_33$condition))

## --- Stage 2: filter + log2 transform, replicating the authors' literal code ---
# Note: this is a much weaker filter than Project 07's convention
# (rowSums(counts >= 10) >= <group size>), which requires a minimum count
# depth per sample within a group, not just a minimum total across all
# samples combined. Matched here to the authors' literal code for
# reproduction fidelity, not because it's the stronger choice -- Project 07's
# filter is the more defensible one for an independent analysis.
keep <- rowSums(counts_raw) > 1
counts_filtered <- counts_raw[keep, ]
cat(sprintf("\nGenes after prefilter (rowSums > 1): %d\n", nrow(counts_filtered)))

sf <- estimateSizeFactorsForMatrix(counts_filtered)
data_scaled_UNUSED <- t(t(counts_filtered) / sf)  # computed, unused — matches authors' code
log2_data <- as.matrix(log2(counts_filtered + 1))

## ALTERNATIVE (presumably-intended, NOT what authors' code literally does):
## log2_data <- as.matrix(log2(data_scaled_UNUSED + 1))

## --- Stage 3: model matrices ---
# Note: mod1 here plays the same role as Project 07's design-formula model
# matrix (~ replicate + co_cultured), but it only has to protect `condition`
# because that's the one known factor -- there's no equivalent of
# `replicate` to block on, since pairing/batch structure in this dataset is
# exactly what SV1 is trying to recover, not something already recorded.
mod1 <- model.matrix(~ condition, data = coldata_33)
mod0 <- model.matrix(~ 1, data = coldata_33)

## --- Stage 4: residualize, then estimate n.sv (matching authors' script) ---
Y.r <- t(resid(lm(t(log2_data) ~ mod1)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
cat(sprintf("\nEstimated number of surrogate variables: %d\n", n.sv))

## --- Stage 5: run SmartSVA with authors' explicit parameters ---
sv_obj <- smartsva.cpp(log2_data, mod = mod1, mod0 = mod0, B = 1000, alpha = 1, n.sv = n.sv)

SV1 <- sv_obj$sv[, 1]
names(SV1) <- colnames(counts_raw)

cat("\nRecomputed SV1 (33 samples):\n")
print(SV1)

## --- Stage 6: validate against authors' published SV1 for the 30 overlapping samples ---
# Note: the design-matrix check done earlier in this project (SV1 negative
# for all PH samples, positive for all POA samples) showed SV1 perfectly
# tracks brain region -- information `condition` already encodes. This is
# the direct analogue of Project 07's write-up prompt 10 ("does a factor's
# name describe what it actually encodes?"): SV1 is nominally a generic
# surrogate variable, but it is partly redundant with condition, the same
# way Project 07 asks you to check whether `replicate` is an honest name
# for what it encodes.
authors_sv1 <- c(
  SRR32154025 = -0.128583969, SRR32154024 = -0.128081602, SRR32154023 = -0.113852846,
  SRR32154022 = -0.125847191, SRR32154021 = -0.166587641, SRR32154020 = -0.165658184,
  SRR32154019 = -0.208081011, SRR32154018 = -0.022100324, SRR32154017 =  0.182329708,
  SRR32154016 =  0.262111161, SRR32154015 =  0.155194487, SRR32154014 =  0.134170043,
  SRR32154013 =  0.207764522, SRR32154012 =  0.179810192, SRR32154011 =  0.056791881,
  SRR32154010 =  0.264846819, SRR32154008 =  0.170742085, SRR32154007 =  0.111604448,
  SRR32154006 =  0.098283414, SRR32154005 =  0.006198624, SRR32154004 =  0.261545312,
  SRR32154003 =  0.025905346, SRR32154002 =  0.227129461, SRR32154001 =  0.016763236,
  SRR32153999 = -0.236061327, SRR32153998 = -0.221466608, SRR32153996 = -0.220708047,
  SRR32153995 = -0.231964587, SRR32153994 = -0.140147975, SRR32153993 = -0.262940444
)
common <- intersect(names(authors_sv1), names(SV1))
correlation <- cor(authors_sv1[common], SV1[common])
cat(sprintf("\nSV1 correlation with authors' GitHub SV1 (%d shared samples): %.4f\n",
            length(common), correlation))

## --- Stage 7: assemble final coldata and save ---
coldata_final <- coldata_33
coldata_final$SV1 <- SV1[rownames(coldata_final)]
coldata_final$condition <- relevel(coldata_final$condition, ref = "PH_intact")

saveRDS(coldata_final, "~/binf6430_repro_env/coldata_final.rds")
saveRDS(counts_filtered, "~/binf6430_repro_env/counts_filtered.rds")

cat("\nSaved coldata_final.rds and counts_filtered.rds — ready for DESeq2.\n")
