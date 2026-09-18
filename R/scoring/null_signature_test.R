#!/usr/bin/env Rscript
# Specificity test: would any six genes do this?
#
# pLSC6's claim is that these six genes in particular track leukemic stemness. The
# sceptical alternative is that any six genes with these weights would move under a
# drug that perturbs the transcriptome. This draws random six-gene signatures, keeps
# the published coefficient vector, and recomputes every contrast, giving an empirical
# null for each observed effect.
#
#   Rscript R/null_signature_test.R [n_permutations]

set.seed(20260917)
source(file.path("R", "lsc_scores.R"))
NPERM <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 1000L }

sig    <- load_signature("pLSC6")
genes  <- sig$gene
coefs  <- sig$coefficient
wanted <- unique(unlist(lapply(seq_len(nrow(sig)), function(i) symbols_for(sig[i, ]))))

ss <- read.csv(file.path("results", "sample_scores.csv"), stringsAsFactors = FALSE)
ss$Line <- ifelse(!is.na(ss$CellLine), ss$CellLine,
            ifelse(!is.na(ss$Genotype), ss$Genotype, ss$GSE))

studies <- list()
for (f in list.files("Toy-Datasets", pattern = "series_matrix\\.txt\\.gz$", full.names = TRUE)) {
  sm <- read_series_matrix(f)
  if (isTRUE(max(sm$expr, na.rm = TRUE) > 50)) sm$expr <- log2(pmax(sm$expr, 0) + 1)
  pm <- try(fetch_platform_annotation(sm$gpl, prefer_symbols = wanted), silent = TRUE)
  if (inherits(pm, "try-error")) next

  # Collapse the WHOLE platform to gene level, so the null can draw from any gene.
  pmk <- pm[!is.na(pm$probe) & !is.na(pm$symbol), ]
  pmk <- pmk[pmk$probe %in% rownames(sm$expr) & nzchar(pmk$symbol) &
             !grepl("///", pmk$symbol, fixed = TRUE), ]
  if (!nrow(pmk)) next
  sub <- sm$expr[pmk$probe, , drop = FALSE]
  mu  <- rowMeans(sub, na.rm = TRUE)
  best <- unlist(tapply(seq_len(nrow(pmk)), pmk$symbol, function(i) i[which.max(mu[i])]))
  M <- sub[best, , drop = FALSE]; rownames(M) <- pmk$symbol[best]
  M <- M[stats::complete.cases(M), , drop = FALSE]

  # Platforms name some signature genes by an alias (GPR56 appears as ADGRG1).
  # Rename those rows to the signature's primary name before anything indexes them.
  for (i in seq_len(nrow(sig))) {
    syn <- symbols_for(sig[i, ])
    hit <- which(rownames(M) %in% syn)
    if (length(hit) && !(sig$gene[i] %in% rownames(M))) rownames(M)[hit[1]] <- sig$gene[i]
  }
  # A platform may not carry every signature gene. Score the observed signature on
  # what is present and give the null the SAME number of genes and the SAME
  # coefficients, so the comparison stays fair.
  have <- genes[genes %in% rownames(M)]
  if (length(have) < 4) { message(sprintf("  %-10s skipped, only %d genes", sm$gse, length(have))); next }
  if (length(have) < length(genes))
    message(sprintf("  %-10s NOTE: %d/%d genes; missing %s", sm$gse, length(have),
                    length(genes), paste(setdiff(genes, have), collapse = ", ")))
  co_s <- setNames(coefs, genes)[have]

  # the null samples only from genes that actually vary
  pool <- which(apply(M, 1, stats::sd) > 0)

  d <- ss[ss$GSE == sm$gse, ]; d <- d[match(colnames(M), d$GSM), ]
  if (all(is.na(d$Agent))) next

  grp <- list()
  for (ln in unique(na.omit(d$Line))) {
    i  <- !is.na(d$Line) & d$Line == ln
    ct <- which(i & !is.na(d$Agent) & d$Agent == "control"); if (length(ct) < 2) next
    for (a in setdiff(unique(na.omit(d$Agent[i])), "control")) {
      tr <- which(i & !is.na(d$Agent) & d$Agent == a); if (length(tr) < 2) next
      grp[[length(grp) + 1]] <- list(key = paste(sm$gse, ln, a, sep = "|"),
                                     agent = a, line = ln, tr = tr, ct = ct)
    }
  }
  if (!length(grp)) next
  studies[[sm$gse]] <- list(M = M, grp = grp, have = have, co = co_s, pool = pool)
  message(sprintf("  %-10s %d genes in pool, %d signature genes, %d contrasts",
                  sm$gse, length(pool), length(have), length(grp)))
}

# effect for one signature (row indices) in one study
eff <- function(M, idx, grp, co) {
  z <- as.numeric(scale(as.numeric(t(M[idx, , drop = FALSE]) %*% co)))
  vapply(grp, function(g) mean(z[g$tr]) - mean(z[g$ct]), numeric(1))
}

obs <- unlist(lapply(studies, function(s) eff(s$M, match(s$have, rownames(s$M)), s$grp, s$co)))
keys <- unlist(lapply(studies, function(s) vapply(s$grp, function(g) g$key, character(1))))
names(obs) <- keys

message(sprintf("\nPermuting %d random six-gene signatures ...", NPERM))
null <- matrix(NA_real_, nrow = NPERM, ncol = length(obs), dimnames = list(NULL, keys))
for (p in seq_len(NPERM)) {
  v <- unlist(lapply(studies, function(s)
        eff(s$M, sample(s$pool, length(s$have)), s$grp, s$co)))
  null[p, ] <- v
}

emp <- vapply(seq_along(obs), function(j)
  (1 + sum(abs(null[, j]) >= abs(obs[j]))) / (NPERM + 1), numeric(1))

parts <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
res <- data.frame(GSE = parts[,1], Line = parts[,2], Agent = parts[,3],
                  observed = round(obs, 3), null_p = round(emp, 3),
                  stringsAsFactors = FALSE)
res <- res[order(res$null_p), ]
cat("\n=== Per contrast: how often does a random six-gene signature match it? ===\n\n")
print(res, row.names = FALSE)

cat(sprintf("\nGlobal: mean |effect| observed %.3f vs random %.3f (p = %.4f)\n",
            mean(abs(obs)), mean(abs(null)),
            (1 + sum(rowMeans(abs(null)) >= mean(abs(obs)))) / (NPERM + 1)))

# The headline claim: does Pinometostat separate MLL-rearranged from wild-type lines
# by more than a random signature would?
MLLR <- c("NOMO-1", "THP-1", "MV4-11", "MOLM-13"); MLLW <- c("HL-60", "U-937", "OCI-AML3")
i_r <- which(res$Agent == "Pinometostat" & res$Line %in% MLLR)
i_w <- which(res$Agent == "Pinometostat" & res$Line %in% MLLW)
if (length(i_r) && length(i_w)) {
  jr <- match(res$GSE[i_r], parts[,1]) # column indices in null
  jr <- which(keys %in% paste(res$GSE[i_r], res$Line[i_r], res$Agent[i_r], sep = "|"))
  jw <- which(keys %in% paste(res$GSE[i_w], res$Line[i_w], res$Agent[i_w], sep = "|"))
  sep_obs  <- mean(obs[jr]) - mean(obs[jw])
  sep_null <- rowMeans(null[, jr, drop = FALSE]) - rowMeans(null[, jw, drop = FALSE])
  cat(sprintf("\nPinometostat, MLL-rearranged minus wild type: %.3f\n", sep_obs))
  cat(sprintf("  random signatures achieving a separation this extreme: p = %.4f\n",
              (1 + sum(sep_null <= sep_obs)) / (NPERM + 1)))
}

saveRDS(list(obs = obs, null = null, keys = keys,
             global_obs = mean(abs(obs)), global_null = rowMeans(abs(null)),
             sep_obs = if (exists("sep_obs")) sep_obs else NA_real_,
             sep_null = if (exists("sep_null")) sep_null else NULL),
        file.path("results", "null_drug_test.rds"))
write.csv(res, file.path("results", "null_signature_test.csv"), row.names = FALSE)
cat("\nWrote results/null_signature_test.csv\n")
