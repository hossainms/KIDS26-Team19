#!/usr/bin/env Rscript
# Corpus forest plot: the contrasts that pass FDR, and whether any are signature-specific.
#
# 95 contrasts will not fit one page, and showing all of them buries the point. This plots
# every contrast surviving Benjamini-Hochberg at q < 0.10, then marks which of those also
# exceed a random six-gene signature. The gap between the two columns is the finding.
#
#   Rscript R/figure_corpus_ranking.R

s <- read.csv(file.path("results","corpus_sample_scores.csv"), stringsAsFactors = FALSE)
s$Line <- ifelse(!is.na(s$CellLine), s$CellLine,
           ifelse(!is.na(s$Genotype), s$Genotype, s$GSE))
s <- s[!is.na(s$Agent), ]
s$Z <- ave(s$Score, s$GSE, FUN = function(v)
  if (stats::sd(v, na.rm = TRUE) > 0) as.numeric(scale(v)) else rep(NA_real_, length(v)))

rows <- list()
for (g in unique(s$GSE)) for (ln in unique(s$Line[s$GSE == g])) {
  d  <- s[s$GSE == g & s$Line == ln, ]
  ct <- d$Z[d$Agent == "control"]; if (length(ct) < 2) next
  for (a in setdiff(unique(d$Agent), "control")) {
    tr <- d$Z[d$Agent == a]; if (length(tr) < 2) next
    tt <- try(t.test(tr, ct), silent = TRUE); if (inherits(tt, "try-error")) next
    rows[[length(rows)+1]] <- data.frame(GSE=g, Line=ln, Agent=a,
      nT=length(tr), nC=length(ct), est=unname(diff(rev(tt$estimate))),
      lo=tt$conf.int[1], hi=tt$conf.int[2], p=tt$p.value, stringsAsFactors=FALSE)
  }
}
r <- do.call(rbind, rows); r$q <- p.adjust(r$p, "BH")
TOTAL <- nrow(r)

nt <- read.csv(file.path("results","corpus_null_test.csv"), stringsAsFactors = FALSE)
r$nullp <- nt$null_p[match(paste(r$GSE,r$Line,r$Agent), paste(nt$GSE,nt$Line,nt$Agent))]
r$specific <- !is.na(r$nullp) & r$nullp < 0.05

r <- r[r$q < 0.10, ]
r <- r[!duplicated(paste(r$Agent, r$Line, round(r$est,3))), ]
r <- r[order(r$est), ]
# Labels go through a base PDF device, so keep them ASCII.
ascii <- function(x) { x <- iconv(x, "UTF-8", "ASCII", sub = ""); trimws(x) }
trunc_at <- function(x, n) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 2), ".."), x)
r$Label <- trunc_at(ascii(r$Agent), 22)
r$LineL <- trunc_at(ascii(r$Line), 18)

INK<-"#16191F"; MUT<-"#6B7480"; GRID<-"#E4E8ED"; SPEC<-"#14532D"; PLAIN<-"#93A0AC"
y <- seq_len(nrow(r))

draw <- function() {
  par(mar=c(0,0,0,0), family="sans", xpd=NA)
  xlo <- -3.2; xhi <- 2.4
  cAg <- xlo-8.4; cLn <- xlo-4.6; cN <- xlo-0.5
  cEs <- xhi+0.5; cQ <- xhi+4.2; cP <- xhi+6.1; cS <- xhi+6.7
  yTop <- -2.6; yAx <- max(y)+1.6; yLab <- max(y)+2.5; yTit <- max(y)+3.5; yLeg <- max(y)+4.8
  plot(NA, xlim=c(cAg,cS), ylim=c(yLeg+2.6, yTop-1.2), axes=FALSE, xlab="", ylab="")

  text(cAg,yTop,"Agent",adj=0,font=2,cex=.95,col=INK)
  text(cLn,yTop,"Cell line",adj=0,font=2,cex=.95,col=INK)
  text(cN, yTop,"n",adj=1,font=2,cex=.95,col=INK)
  text(cEs,yTop,expression(bold(paste(Delta," (95% CI)"))),adj=0,cex=.95,col=INK)
  text(cQ, yTop,"BH q",adj=1,font=2,cex=.95,col=INK)
  text(cP, yTop,"null P",adj=1,font=2,cex=.95,col=INK)
  segments(cAg,yTop+.7,cS,yTop+.7,col=INK,lwd=.9)
  text(-1.6,yTop+1.8,"lower stemness",cex=.86,col=MUT,adj=.5)
  text( 1.4,yTop+1.8,"higher",cex=.86,col=MUT,adj=.5)

  at <- seq(-3,2,1)
  segments(at,yTop+2.4,at,yAx,col=GRID,lwd=.6); segments(0,yTop+2.4,0,yAx,col=INK,lwd=.9)

  for (i in seq_len(nrow(r))) {
    col <- if (r$specific[i]) SPEC else PLAIN
    x1 <- max(r$lo[i],xlo); x2 <- min(r$hi[i],xhi)
    segments(x1,y[i],x2,y[i],col=col,lwd=1.1)
    if (r$lo[i]<xlo) arrows(x1,y[i],xlo-.14,y[i],length=.03,angle=22,col=col,lwd=1.1)
    else segments(x1,y[i]-.2,x1,y[i]+.2,col=col,lwd=1.1)
    if (r$hi[i]>xhi) arrows(x2,y[i],xhi+.14,y[i],length=.03,angle=22,col=col,lwd=1.1)
    else segments(x2,y[i]-.2,x2,y[i]+.2,col=col,lwd=1.1)
    points(r$est[i],y[i],pch=if(r$specific[i]) 22 else 0,bg=col,col=col,cex=1.05,lwd=1.1)
    text(cAg,y[i],r$Label[i],adj=0,cex=.88,col=INK,font=if(r$specific[i]) 2 else 1)
    text(cLn,y[i],r$LineL[i],adj=0,cex=.86,col=INK)
    text(cN, y[i],paste0(r$nT[i],"/",r$nC[i]),adj=1,cex=.84,col=MUT)
    text(cEs,y[i],sprintf("%+.2f (%.2f, %.2f)",r$est[i],r$lo[i],r$hi[i]),adj=0,cex=.84,col=INK)
    text(cQ, y[i],sprintf("%.3f",r$q[i]),adj=1,cex=.84,col=INK)
    text(cP, y[i],if(is.na(r$nullp[i])) "-" else sprintf("%.3f",r$nullp[i]),
         adj=1,cex=.84,col=if(r$specific[i]) SPEC else MUT,font=if(r$specific[i]) 2 else 1)
    if (r$specific[i]) points(cS,y[i],pch=16,cex=.7,col=SPEC)
  }

  segments(xlo,yAx,xhi,yAx,col=INK,lwd=.8); segments(at,yAx,at,yAx+.22,col=INK,lwd=.8)
  text(at,yLab,at,cex=.88,col=INK)
  text(mean(c(xlo,xhi)),yTit,"Standardised shift in pLSC6 versus control (within-study SD)",
       cex=1.0,col=INK)
  segments(cAg,max(y)+1.0,cS,max(y)+1.0,col=INK,lwd=.9)

  points(cAg,yLeg,pch=22,bg=SPEC,col=SPEC,cex=1.2)
  text(cAg+.25,yLeg,"exceeds a random six-gene signature",adj=0,cex=.86,col=INK)
  points(cAg+7.6,yLeg,pch=0,col=PLAIN,cex=1.2)
  text(cAg+7.85,yLeg,"does not",adj=0,cex=.86,col=INK)

  fn <- c(
    sprintf("%d of %d contrasts survive Benjamini-Hochberg at q < 0.10 on effect size; all %d are shown.",
            nrow(r), TOTAL, nrow(r)),
    sprintf("Filled points exceed a random six-gene signature at raw P < 0.05. Only %d of %d contrasts do, which is the rate expected by chance,",
            sum(nt$null_p < .05), TOTAL),
    "and after correcting the permutation P values only Pinometostat survives (q = 0.095). Globally pLSC6 does not beat random sets (P = 0.99).",
    "Corpus: 38 GEO series, 613 samples, 41 agents. Many arms are n = 2 and half the studies score on 5 of 6 genes.")
  for (k in seq_along(fn)) text(cAg, yLeg+1.0+(k-1)*0.85, fn[k], adj=0, cex=.8, col=MUT)
}

H <- max(5.2, 1.6 + nrow(r)*0.30)
pdf(file.path("results","figure_corpus_ranking.pdf"), width=12.4, height=H, useDingbats=FALSE)
draw(); invisible(dev.off())
png(file.path("results","figure_corpus_ranking.png"), width=12.4, height=H, units="in", res=600)
draw(); invisible(dev.off())
cat(sprintf("%d of %d contrasts plotted (q<0.10); %d signature-specific\n",
            nrow(r), TOTAL, sum(r$specific)))
