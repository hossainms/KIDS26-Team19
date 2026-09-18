#!/usr/bin/env Rscript
# Where the signature carries specific information, and where it does not.
#
# Three permutation tests, each against random six-gene signatures carrying pLSC6's
# published coefficients. The contrast between panel A and panel B is the result: the
# signature is specifically prognostic in patients and largely non-specific as a
# pharmacodynamic readout in cell lines.

suppressPackageStartupMessages(library(survival))

sv <- readRDS(file.path("results", "null_survival_test.rds"))
dg <- readRDS(file.path("results", "corpus_null_draws.rds"))

INK <- "#16191F"; MUT <- "#6B7480"; GRID <- "#E4E8ED"
NULLC <- "#C7D0DA"; NULLE <- "#9AA6B2"; HIT <- "#14532D"; MISS <- "#A8202B"

panel <- function(null, obs, side, xlab, title, lab, pval, hit) {
  h <- hist(null, breaks = 34, plot = FALSE)
  xr <- range(c(h$breaks, obs)) + c(-.04, .04) * diff(range(c(h$breaks, obs)))
  plot(NA, xlim = xr, ylim = c(0, max(h$counts) * 1.30), axes = FALSE, xlab = "", ylab = "")

  tail_sel <- if (side == "right") h$mids >= obs else h$mids <= obs
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts,
       col = ifelse(tail_sel, if (hit) HIT else MISS, NULLC),
       border = NULLE, lwd = .4)

  axis(1, col = INK, col.axis = INK, cex.axis = .95, lwd = .8, tck = -.022, mgp = c(3, .7, 0))
  axis(2, col = INK, col.axis = INK, cex.axis = .95, lwd = .8, las = 1,
       tck = -.022, mgp = c(3, .6, 0))
  mtext(xlab, side = 1, line = 2.5, cex = .82, col = INK)
  mtext("random signatures", side = 2, line = 2.9, cex = .82, col = INK)

  segments(obs, 0, obs, max(h$counts) * 1.10, col = if (hit) HIT else MISS, lwd = 2.4)
  points(obs, max(h$counts) * 1.10, pch = 25, bg = if (hit) HIT else MISS,
         col = if (hit) HIT else MISS, cex = 1.15)
  text(obs, max(h$counts) * 1.20, lab, col = if (hit) HIT else MISS,
       cex = .95, font = 2,
       adj = if (obs > mean(xr)) 1.03 else -0.03)
  text(obs, max(h$counts) * 1.20, sprintf("  P = %s  ", pval), col = INK, cex = .88,
       adj = if (obs > mean(xr)) c(1.03, 1.9) else c(-0.03, 1.9))

  mtext(title, side = 3, line = 1.25, adj = 0, cex = .98, font = 2, col = INK)
}

draw <- function() {
  par(mfrow = c(1, 3), mar = c(4.4, 4.6, 4.2, 1.4), family = "sans",
      oma = c(0.4, 0.4, 2.6, 0.4))

  panel(sv$null[, "z"], sv$obs["z"], "right",
        "Cox z statistic", "A   Survival, 351 TARGET patients",
        "pLSC6 4.43", "0.005", TRUE)
  mtext("31% of random signatures reach P < 0.05 here", side = 3,
        line = 0.1, adj = 0, cex = .72, col = MUT)

  panel(dg$global_null, dg$global_obs, "right",
        "mean absolute effect, 95 contrasts", "B   Drug response, 38 GEO series",
        "pLSC6 0.58", "0.99", FALSE)
  mtext("95 contrasts, 41 agents: random six-gene sets do better", side = 3,
        line = 0.1, adj = 0, cex = .72, col = MUT)

  panel(dg$sep_null, dg$sep_obs, "left",
        "MLL-rearranged minus wild type", "C   Pinometostat, genotype split",
        "-1.60", "0.004", TRUE)
  mtext("the one cell-line result that is signature-specific", side = 3,
        line = 0.1, adj = 0, cex = .72, col = MUT)

  mtext("pLSC6 against 1,000 random six-gene signatures carrying the same coefficients",
        side = 3, outer = TRUE, line = 0.5, adj = 0, cex = .92, col = INK, font = 2)
}

pdf(file.path("results", "figure_null_tests.pdf"), width = 12.6, height = 4.5,
    useDingbats = FALSE); draw(); invisible(dev.off())
png(file.path("results", "figure_null_tests.png"), width = 12.6, height = 4.5,
    units = "in", res = 600); draw(); invisible(dev.off())
cat("Wrote results/figure_null_tests.{pdf,png}\n")
