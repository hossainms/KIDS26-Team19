#!/usr/bin/env Rscript
# Every contrast that exceeded a random six-gene signature, across both corpora.
#
# 257 contrasts will not fit a page, and plotting only the ones significant on effect size
# would show 80-odd rows that say nothing about specificity. This plots the 18 that beat
# their permutation null at raw P < 0.05 -- against a chance expectation of 12.9 -- so the
# reader can see both how few there are and what they consist of.

ar <- read.csv(file.path("results","corpus_drug_ranking.csv"),  stringsAsFactors = FALSE)
an <- read.csv(file.path("results","corpus_null_test.csv"),     stringsAsFactors = FALSE)
nr <- read.csv(file.path("results","ncbi_drug_ranking.csv"),    stringsAsFactors = FALSE)
nn <- read.csv(file.path("results","ncbi_null_test.csv"),       stringsAsFactors = FALSE)
key <- function(d) paste(d$GSE, d$Line, d$Agent)
ar$src <- "array"; nr$src <- "counts"
ar$nullp <- an$null_p[match(key(ar), key(an))]
nr$nullp <- nn$null_p[match(key(nr), key(nn))]
cols <- c("GSE","Line","Agent","nT","nC","est","p","q","nullp","src")
r <- rbind(ar[,cols], nr[,cols])
TOTAL <- nrow(r); EXP <- 0.05 * TOTAL
r$qall <- p.adjust(r$p, "BH")
r <- r[!is.na(r$nullp) & r$nullp < 0.05, ]
r <- r[order(r$est), ]

ascii <- function(x) trimws(iconv(x, "UTF-8", "ASCII", sub = ""))
trunc_at <- function(x, n) ifelse(nchar(x) > n, paste0(substr(x, 1, n-2), ".."), x)
r$Label <- trunc_at(ascii(r$Agent), 26); r$LineL <- trunc_at(ascii(r$Line), 16)
ctrl <- grepl("scramble|shRNA control|vehicle|DMSO", r$Agent, ignore.case = TRUE)

INK<-"#16191F"; MUT<-"#6B7480"; GRID<-"#E4E8ED"; ARR<-"#14532D"; CNT<-"#1E5A8A"; FLAG<-"#A8202B"
y <- seq_len(nrow(r))

draw <- function() {
  par(mar = c(0,0,0,0), family = "sans", xpd = NA)
  xlo <- -3.0; xhi <- 2.8
  cAg <- xlo-9.0; cLn <- xlo-5.0; cN <- xlo-1.6; cSr <- xlo-0.5
  cEs <- xhi+0.5; cQ <- xhi+3.1; cP <- xhi+4.8
  yTop <- -2.6; yAx <- max(y)+1.6; yLab <- max(y)+2.5; yTit <- max(y)+3.5; yLeg <- max(y)+4.9
  plot(NA, xlim = c(cAg,cP), ylim = c(yLeg+3.0, yTop-1.2), axes = FALSE, xlab="", ylab="")

  text(cAg,yTop,"Agent",adj=0,font=2,cex=.95,col=INK)
  text(cLn,yTop,"Cell line",adj=0,font=2,cex=.95,col=INK)
  text(cN, yTop,"n",adj=1,font=2,cex=.95,col=INK)
  text(cSr,yTop,"source",adj=1,font=2,cex=.95,col=INK)
  text(cEs,yTop,expression(bold(paste(Delta," (SD)"))),adj=0,cex=.95,col=INK)
  text(cQ, yTop,"BH q",adj=1,font=2,cex=.95,col=INK)
  text(cP, yTop,"null P",adj=1,font=2,cex=.95,col=INK)
  segments(cAg,yTop+.7,cP,yTop+.7,col=INK,lwd=.9)
  at <- seq(-3,2,1)
  segments(at,yTop+1.6,at,yAx,col=GRID,lwd=.6); segments(0,yTop+1.6,0,yAx,col=INK,lwd=.9)

  for (i in seq_len(nrow(r))) {
    col <- if (ctrl[i]) FLAG else if (r$src[i]=="array") ARR else CNT
    segments(0,y[i],r$est[i],y[i],col=col,lwd=1.2)
    points(r$est[i],y[i],pch=22,bg=col,col=col,cex=1.0,lwd=1.1)
    text(cAg,y[i],r$Label[i],adj=0,cex=.88,col=if (ctrl[i]) FLAG else INK,
         font=if (ctrl[i]) 2 else 1)
    text(cLn,y[i],r$LineL[i],adj=0,cex=.86,col=INK)
    text(cN, y[i],paste0(r$nT[i],"/",r$nC[i]),adj=1,cex=.84,col=MUT)
    text(cSr,y[i],r$src[i],adj=1,cex=.82,col=MUT)
    text(cEs,y[i],sprintf("%+.2f",r$est[i]),adj=0,cex=.86,col=INK)
    text(cQ, y[i],sprintf("%.3f",r$qall[i]),adj=1,cex=.84,col=MUT)
    text(cP, y[i],sprintf("%.3f",r$nullp[i]),adj=1,cex=.84,col=INK)
  }

  segments(xlo,yAx,xhi,yAx,col=INK,lwd=.8); segments(at,yAx,at,yAx+.22,col=INK,lwd=.8)
  text(at,yLab,at,cex=.88,col=INK)
  text(mean(c(xlo,xhi)),yTit,"Standardised shift in pLSC6 versus control (within-study SD)",cex=1.0,col=INK)
  segments(cAg,max(y)+1.0,cP,max(y)+1.0,col=INK,lwd=.9)

  points(cAg,yLeg,pch=22,bg=ARR,col=ARR,cex=1.2); text(cAg+.25,yLeg,"array series matrix",adj=0,cex=.86,col=INK)
  points(cAg+5.0,yLeg,pch=22,bg=CNT,col=CNT,cex=1.2); text(cAg+5.25,yLeg,"NCBI RNA-seq counts",adj=0,cex=.86,col=INK)
  points(cAg+10.6,yLeg,pch=22,bg=FLAG,col=FLAG,cex=1.2); text(cAg+10.85,yLeg,"not a drug",adj=0,cex=.86,col=INK)

  fn <- c(sprintf("All %d contrasts (of %d across 76 series) that exceed a random six-gene signature at raw permutation P < 0.05.",
                  nrow(r), TOTAL),
          sprintf("Chance alone predicts %.1f. None survives correction of the permutation P values, and neither corpus beats random overall (P = 0.99, P = 0.39).", EXP),
          "The highest-ranking contrast is a scramble shRNA negative control, which indicates what this tail is made of.")
  for (k in seq_along(fn)) text(cAg, yLeg+1.1+(k-1)*0.85, fn[k], adj=0, cex=.8, col=MUT)
}

H <- max(5.0, 2.2 + nrow(r)*0.34)
pdf(file.path("results","figure_combined_ranking.pdf"), width=13.0, height=H, useDingbats=FALSE)
draw(); invisible(dev.off())
png(file.path("results","figure_combined_ranking.png"), width=13.0, height=H, units="in", res=600)
draw(); invisible(dev.off())
cat(sprintf("%d contrasts beat the null of %d total (chance %.1f); %d are negative controls\n",
            nrow(r), TOTAL, EXP, sum(ctrl)))
