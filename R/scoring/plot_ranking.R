#!/usr/bin/env Rscript
# Draw the agent ranking split by cell line.
#
# The collapsed view averages a drug's effect across every line it was tested in, which
# can cancel two large opposite effects to nothing. Splitting restores them, and shows
# that the largest responses track whether the line carries the drug's target dependency.

rk <- read.csv(file.path("results", "drug_ranking.csv"), stringsAsFactors = FALSE)
rk$CellLine[is.na(rk$CellLine)] <- rk$GSE[is.na(rk$CellLine)]

# Order agents by their strongest single response, and lines within an agent by effect.
best  <- tapply(rk$DeltaZ, rk$Agent, min)
rk$ag <- factor(rk$Agent, levels = names(sort(best, decreasing = TRUE)))
rk    <- rk[order(rk$ag, -rk$DeltaZ), ]

# On target = the line carries the dependency the drug is built against.
on_target <- grepl("rearranged|NfsB-positive", rk$Context, ignore.case = TRUE)
off_target <- grepl("wild", rk$Context, ignore.case = TRUE)

ON  <- "#1B5E44"; OFF <- "#9AA6B2"; UNK <- "#4E8C73"; POS <- "#A8202B"
col <- ifelse(rk$DeltaZ > 0 & off_target, POS,
       ifelse(on_target, ON, ifelse(off_target, OFF, UNK)))
col[rk$DeltaZ > 0 & !off_target] <- POS

# One slot per contrast, a gap between agents.
grp <- as.integer(rk$ag)
y   <- seq_len(nrow(rk)) + cumsum(c(0, diff(grp) != 0)) * 0.9

png(file.path("results", "archive", "drug_ranking_by_cellline.png"),
    width = 1900, height = 1250, res = 165)
op <- par(mar = c(5.2, 11.5, 5.6, 8.5), xpd = NA)

xr <- range(c(rk$DeltaZ, 0)) + c(-0.45, 0.35)
plot(NA, xlim = xr, ylim = c(max(y) + 1.2, min(y) - 2.6), axes = FALSE,
     xlab = "", ylab = "")

abline(v = pretty(xr), col = "#E8ECF1", lwd = 1)
rect(0, y - 0.38, rk$DeltaZ, y + 0.38, col = col, border = NA)
abline(v = 0, col = "#3C4654", lwd = 1.6)

axis(1, at = pretty(xr), col = "#B8C1CC", col.axis = "#3C4654", cex.axis = 0.88)
mtext("shift in pLSC6 vs the control arm of the same study  (within-study SD)",
      side = 1, line = 2.9, cex = 0.95, col = "#3C4654")

# cell line labels
text(xr[1] + 0.06, y, rk$CellLine, adj = 0, cex = 0.82, col = "#20262F")
# n and context on the right
text(xr[2], y, sprintf("n=%d/%d", rk$nTreated, rk$nControl),
     adj = 1, cex = 0.72, col = "#8892A2")

# agent group labels, centred on each block
for (a in levels(rk$ag)) {
  i <- which(rk$ag == a)
  yy <- mean(y[i])
  text(xr[1] - 0.30, yy, a, adj = 1, cex = 0.95, font = 2, col = "#20262F")
  if (length(i) > 1)
    lines(rep(xr[1] - 0.22, 2), c(min(y[i]) - 0.42, max(y[i]) + 0.42),
          col = "#C3CBD8", lwd = 2)
}

title(main = "Every drug, every cell line", adj = 0, line = 3.5,
      cex.main = 1.45, col.main = "#12161C", font.main = 2)
mtext("20 treated-vs-control contrasts across 4 studies. Left of the line means the drug moved cells away from the stem-like state.",
      side = 3, line = 2.1, adj = 0, cex = 0.86, col = "#5A6475")

legend("top", horiz = TRUE, bty = "n", inset = c(0, -0.085), cex = 0.8,
       fill = c(ON, UNK, OFF, POS), border = NA,
       legend = c("carries the drug's target dependency", "dependency not recorded",
                  "target-negative line", "score increased"))
par(op); dev.off()
cat("Wrote results/drug_ranking_by_cellline.png\n")

# the contrast the collapsed chart destroys
cat("\nWhat splitting recovers:\n")
for (a in c("Pinometostat", "LSD1 OG86")) {
  s <- rk[rk$Agent == a, ]
  if (!nrow(s)) next
  on <- grepl("rearranged|NfsB-positive", s$Context, ignore.case = TRUE)
  cat(sprintf("  %-13s target-positive %+.2f   target-negative %+.2f   collapsed %+.2f\n",
              a, mean(s$DeltaZ[on]), mean(s$DeltaZ[!on]), mean(s$DeltaZ)))
}
