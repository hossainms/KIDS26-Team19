#!/usr/bin/env Rscript
# SUPERSEDED by R/figure_corpus_ranking.R, which runs on 95 corpus contrasts rather than
# the 20 pilot contrasts here. Output goes to results/archive/. Kept for provenance.
#
# Publication figure: forest plot of pLSC6 response by agent and cell line.
#
# Each row is one treated-versus-control contrast within a single study. Effects are
# standardised within study, so they are comparable across platforms. Intervals are
# Welch t 95% confidence limits; most arms are n = 2, so the intervals are wide and
# are the point of the figure rather than a footnote.
#
#   Rscript R/figure_ranking.R
# Writes results/figure_drug_ranking.{pdf,png}

s <- read.csv(file.path("results", "sample_scores.csv"), stringsAsFactors = FALSE)
s$Line <- ifelse(!is.na(s$CellLine), s$CellLine,
           ifelse(!is.na(s$Genotype), s$Genotype, s$GSE))
s <- s[!is.na(s$Agent), ]
s$Z <- ave(s$Score, s$GSE, FUN = function(v) as.numeric(scale(v)))

rows <- list()
for (g in unique(s$GSE)) for (ln in unique(s$Line[s$GSE == g])) {
  d  <- s[s$GSE == g & s$Line == ln, ]
  ct <- d$Z[d$Agent == "control"]
  if (length(ct) < 2) next
  for (a in setdiff(unique(d$Agent), "control")) {
    tr <- d$Z[d$Agent == a]
    if (length(tr) < 2) next
    tt  <- t.test(tr, ct)
    ctx <- unique(na.omit(d$Genotype))
    rows[[length(rows) + 1]] <- data.frame(
      GSE = g, Line = ln, Agent = a, nT = length(tr), nC = length(ct),
      est = unname(diff(rev(tt$estimate))), lo = tt$conf.int[1], hi = tt$conf.int[2],
      p = tt$p.value,
      Context = if (length(ctx)) ctx[1]
                else if (grepl("NfsB", ln, ignore.case = TRUE)) "target-positive"
                else if (grepl("wild", ln, ignore.case = TRUE)) "target-negative"
                else NA_character_,
      stringsAsFactors = FALSE)
  }
}
r <- do.call(rbind, rows)

r$q <- p.adjust(r$p, "BH")

# Contrasts that beat a random six-gene signature (R/null_signature_test.R).
nf <- file.path("results", "null_signature_test.csv")
r$nullp <- NA_real_
if (file.exists(nf)) {
  nt <- read.csv(nf, stringsAsFactors = FALSE)
  key  <- paste(r$GSE, r$Line, r$Agent);  nkey <- paste(nt$GSE, nt$Line, nt$Agent)
  r$nullp <- nt$null_p[match(key, nkey)]
}
r$specific <- !is.na(r$nullp) & r$nullp < 0.05

# Two platforms carry no FAM30A, so those studies score on 5 of the 6 genes.
r$partial <- r$GSE %in% c("GSE144638", "GSE155640")

r$onTarget <- grepl("rearranged|target-positive", r$Context, ignore.case = TRUE)
r$offTarget <- grepl("wild|target-negative", r$Context, ignore.case = TRUE)
r$sig <- r$p < 0.05

# Display labels only. The full strings stay in contrast_statistics.csv; the figure
# uses the compact genotype notation so the column cannot run into the n column.
r$Label <- r$Line
r$Label <- sub("^NfsB-positive THP1 cells$", "THP1 NfsB+", r$Label)
r$Label <- sub("^THP1 wild-type cells$",     "THP1 WT",    r$Label)
r$Label <- sub("\\s+cells$", "", r$Label)

best <- tapply(r$est, r$Agent, min)
r$ag <- factor(r$Agent, levels = names(sort(best, decreasing = TRUE)))
r <- r[order(r$ag, -r$est), ]

grp <- as.integer(r$ag)
y   <- (seq_len(nrow(r)) + cumsum(c(0, diff(grp) != 0)) * 1.15) * 1.18

INK <- "#16191F"; MUT <- "#6B7480"; GRID <- "#E4E8ED"
ON  <- "#14532D"; OFFC <- "#8B95A1"; NEU <- "#2F6F52"

draw <- function() {
  par(mar = c(0, 0, 0, 0), family = "sans", xpd = NA)

  # Plot geometry in data units. Columns are placed explicitly so nothing collides.
  xlo <- -6.3; xhi <- 1.9
  cAg <- xlo - 8.0      # agent
  cLn <- xlo - 5.0      # cell line
  cN  <- xlo - 0.6      # n, right aligned
  cEs <- xhi + 0.6      # estimate (95% CI)
  cP  <- xhi + 4.6      # P, right aligned
  cQ  <- xhi + 6.2      # BH q, right aligned
  cS  <- xhi + 6.9      # specificity mark

  yTop <- -4.4                       # header row
  yAx  <- max(y) + 1.7               # axis line
  yLab <- max(y) + 2.7               # tick labels
  yTit <- max(y) + 3.9               # axis title
  yLeg <- max(y) + 5.3               # legend

  plot(NA, xlim = c(cAg, cS), ylim = c(yLeg + 4.3, yTop - 1.4),
       axes = FALSE, xlab = "", ylab = "")

  # header
  text(cAg, yTop, "Agent",     adj = 0, font = 2, cex = .95, col = INK)
  text(cLn, yTop, "Cell line", adj = 0, font = 2, cex = .95, col = INK)
  text(cN,  yTop, "n",         adj = 1, font = 2, cex = .95, col = INK)
  text(cEs, yTop, expression(bold(paste(Delta, " (95% CI)"))), adj = 0, cex = .95, col = INK)
  text(cP,  yTop, "P",   adj = 1, font = 2, cex = .95, col = INK)
  text(cQ,  yTop, "q",   adj = 1, font = 2, cex = .95, col = INK)
  segments(cAg, yTop + 0.75, cS, yTop + 0.75, col = INK, lwd = .9)

  # direction cues, inside the plotting band
  text(-3.2, yTop + 1.9, "lower stemness", cex = .86, col = MUT, adj = .5)
  text( 0.95, yTop + 1.9, "higher",        cex = .86, col = MUT, adj = .5)

  # grid and zero line
  at <- seq(-6, 1, 1)
  segments(at, yTop + 2.6, at, yAx, col = GRID, lwd = .6)
  segments(0,  yTop + 2.6, 0,  yAx, col = INK,  lwd = .9)

  for (i in seq_len(nrow(r))) {
    col <- if (r$onTarget[i]) ON else if (r$offTarget[i]) OFFC else NEU
    x1 <- max(r$lo[i], xlo); x2 <- min(r$hi[i], xhi)
    segments(x1, y[i], x2, y[i], col = col, lwd = 1.1)
    if (r$lo[i] < xlo) arrows(x1, y[i], xlo - .18, y[i], length = .03, angle = 22,
                              col = col, lwd = 1.1)
    else segments(x1, y[i] - .2, x1, y[i] + .2, col = col, lwd = 1.1)
    if (r$hi[i] > xhi) arrows(x2, y[i], xhi + .18, y[i], length = .03, angle = 22,
                              col = col, lwd = 1.1)
    else segments(x2, y[i] - .2, x2, y[i] + .2, col = col, lwd = 1.1)

    se <- (r$hi[i] - r$lo[i]) / 3.92
    points(r$est[i], y[i], pch = if (r$sig[i]) 22 else 0, bg = col, col = col,
           cex = min(0.6 + 0.35 / max(se, 0.08), 1.45), lwd = 1.0)

    text(cLn, y[i], r$Label[i], adj = 0, cex = .90, col = INK)
    text(cN,  y[i], paste0(r$nT[i], "/", r$nC[i]), adj = 1, cex = .86, col = MUT)
    text(cEs, y[i], sprintf("%+.2f (%.2f, %.2f)", r$est[i], r$lo[i], r$hi[i]),
         adj = 0, cex = .86, col = INK)
    text(cP,  y[i], if (r$p[i] < .001) "<0.001" else sprintf("%.3f", r$p[i]),
         adj = 1, cex = .86, col = if (r$sig[i]) INK else MUT)
    text(cQ,  y[i], sprintf("%.2f", r$q[i]), adj = 1, cex = .86,
         col = if (r$q[i] < .10) INK else MUT)
    if (r$specific[i]) text(cS, y[i], "\u25CF", adj = 1, cex = .8, col = ON)
    if (r$partial[i])  text(cLn - 0.28, y[i], "\u2020", adj = 1, cex = .8, col = MUT)
  }

  for (a in levels(r$ag)) {
    i <- which(r$ag == a)
    text(cAg, mean(y[i]), a, adj = 0, cex = .98, font = 2, col = INK)
  }

  # axis
  segments(xlo, yAx, xhi, yAx, col = INK, lwd = .8)
  segments(at, yAx, at, yAx + 0.26, col = INK, lwd = .8)
  text(at, yLab, at, cex = .88, col = INK)
  text(mean(c(xlo, xhi)), yTit,
       "Standardised shift in pLSC6 versus control (within-study SD)",
       cex = 1.02, col = INK)

  # legend, on its own row below the axis title
  lx <- c(cAg, cAg + 5.2, cAg + 10.2)
  points(lx[1], yLeg, pch = 22, bg = ON, col = ON, cex = 1.3)
  text(lx[1] + 0.2, yLeg, "target-dependent line", adj = 0, cex = .86, col = INK)
  points(lx[2], yLeg, pch = 22, bg = OFFC, col = OFFC, cex = 1.3)
  text(lx[2] + 0.2, yLeg, "target-negative line", adj = 0, cex = .86, col = INK)
  points(lx[3], yLeg, pch = 0, col = NEU, cex = 1.3)
  text(lx[3] + 0.2, yLeg, "open symbol: raw P \u2265 0.05", adj = 0, cex = .86, col = INK)

  fn <- c(
    "q, Benjamini-Hochberg over 20 contrasts; no contrast reaches q < 0.05.",
    "\u25CF marks the two contrasts exceeding a random six-gene signature (P < 0.05, 1,000 permutations).",
    "\u2020 GPL23159 and GPL17586 carry no FAM30A, so these studies score on 5 of the 6 genes.",
    "In GSE134742 the agents were given at 0.4, 1.0 and 10 \u00B5M, so cross-agent potency is not comparable.")
  for (k in seq_along(fn))
    text(cAg, yLeg + 0.95 + (k - 1) * 0.78, fn[k], adj = 0, cex = .76, col = MUT)
}

dir.create("results", showWarnings = FALSE)
pdf(file.path("results", "archive", "figure_drug_ranking.pdf"), width = 11.6, height = 8.4,
    useDingbats = FALSE); draw(); invisible(dev.off())
png(file.path("results", "archive", "figure_drug_ranking.png"), width = 11.6, height = 8.4,
    units = "in", res = 600); draw(); invisible(dev.off())

write.csv(r[, c("GSE","Line","Agent","Context","nT","nC","est","lo","hi","p")],
          file.path("results", "contrast_statistics.csv"), row.names = FALSE)
cat(sprintf("%d contrasts, %d significant at P<0.05\n", nrow(r), sum(r$sig)))
cat("Wrote results/figure_drug_ranking.{pdf,png} and contrast_statistics.csv\n")
