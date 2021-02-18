args <- list(vcf.file = '/Users/leraman/PhD/BAF/batch1-gatk-haplotype-joint-annotated.vcf.gz',
             
             singles = c('D1810466', 'D1810468', 'D1810470'), # use for separate sample testing; will detect LOHs
             trio = c(), # use for trio analysis: father, mother, child; will detect LOH and increased MEs
             
             output.dir = '/Users/leraman/PhD/BAF/output2',
             
             cairo.bitmap = F, # some systems require cairo bitmap for plotting
             
             min.seq.depth = 15, # consider variants with at least a depth of min.seq.depth
             
             hist.bins = 100, # number of bins in beta value histograms
             
             bin.size = 5, # in MB, for LOH/ME calling
             min.support = 1, # number of variants required in bin for P-value calculation; for LOH/ME calling
             successive.bins.required = 3, # successive 'significant' bins required for actual significance of middle bin; for LOH/ME calling; odd number
             
             cytoband = '/Users/leraman/PhD/BAF/cyotband.hg38.txt' # if given, filters out gvar and stalk
)

if (length(args$singles) & length(args$trio)){
  stop('Provide single unrelated samples OR a trio of samples\n')
}

if (length(args$trio) & length(args$trio) != 3){
  stop('A trio requires exactly three samples\n')
}

###
# initialize
###

suppressMessages(library('vcfR')) # v1.10.0
suppressMessages(library('GenomicRanges')) # v1.34.0
suppressMessages(library('trio')) # v3.20.0

options(scipen=999)
options(digits=10)
if (args$cairo) options(bitmaptype='cairo')

chr.lengths <- list('chr1'=248956422, 'chr10'=133797422, 'chr11'=135086622, 'chr12'=133275309, 'chr13'=114364328, 'chr14'=107043718,
                    'chr15'=101991189, 'chr16'=90338345, 'chr17'=83257441, 'chr18'=80373285, 'chr19'=58617616, 'chr2'=242193529,
                    'chr20'=64444167, 'chr21'=46709983, 'chr22'=50818468, 'chr3'=198295559, 'chr4'=190214555, 'chr5'=181538259,
                    'chr6'=170805979, 'chr7'=159345973, 'chr8'=145138636, 'chr9'=138394717, 'chrX'=156040895, 'chrY'=57227415) # hg38

chr.order = paste0('chr', c(1:22, 'X'))
chr.lengths <- chr.lengths[chr.order]
chr.cumsum <- c(1, cumsum(chr.lengths[-length(chr.lengths)]) + 1)
names(chr.cumsum) <- chr.order
chr.mids <- chr.cumsum + unlist(chr.lengths) / 2

color.A = rgb(0, 122, 104, maxColorValue = 255)
color.B = rgb(226, 191, 121, maxColorValue = 255)
color.C = rgb(150, 80, 34, maxColorValue = 255)

###
# load
###

load.samples <- function(args){ # https://github.com/leraman/hopla/blob/master/hopla.R
  cat('Loading vcf.gz ...\n')
  vcf <- read.vcfR(args$vcf.file, verbose = F)
  
  ## vcf.A (annot)
  vcf.A <- as.data.frame(vcf@fix)
  if (substr(vcf.A$CHROM[1], 1, 3) != 'chr') vcf.A$CHROM <- paste0('chr', vcf.A$CHROM)
  vcf.A$ID <- as.character(paste0('id', 1:nrow(vcf.A)))
  
  vcf.A <- vcf.A[,c('CHROM', 'POS', 'ID', 'REF', 'ALT')]
  vcf.A$POS <- as.numeric(as.character(vcf.A$POS))
  for (x in c('CHROM', 'ID', 'REF', 'ALT')) vcf.A[[x]] <- as.character(vcf.A[[x]])
  
  hard.mask <- nchar(vcf.A$REF) == 1 & nchar(vcf.A$ALT) == 1 & vcf.A$CHROM %in% chr.order
  if (length(args$cytoband)){
    cyto <- read.csv(args$cytoband, sep = '\t', header = F, stringsAsFactors = F)
    cyto <- cyto[cyto$V5 %in% c('acen', 'gvar', 'stalk'),]
    for (i in 1:nrow(cyto)){
      cyto.mask = vcf.A$CHROM == cyto$V1[i] & vcf.A$POS > cyto$V2[i] - 1e6 & vcf.A$POS < cyto$V3[i] + 1e6
      hard.mask[cyto.mask] <- F
    }
  }
  
  ## vcf.B (data)
  cat('Parsing variants, working ...\n')
  vcf.B.tmp <- as.data.frame(vcf@gt)
  cnames <- strsplit(as.character(vcf.B.tmp$FORMAT[1]), ':')[[1]]
  vcfs <- list()
  once = T
  for (sample in c(args$singles, args$trio)){
    cat(paste0('  ... at ', sample, '\n'))
    vcf.B <- suppressWarnings(as.data.frame(do.call(rbind, strsplit(as.character(vcf.B.tmp[[sample]]), ':')))[,1:length(cnames)])
    colnames(vcf.B) <- cnames
    
    vcf.B <- vcf.B[,c('GT', 'AD', 'DP')]
    vcf.B$DP <- suppressWarnings(as.numeric(as.character(vcf.B$DP)))
    for (x in c('GT', 'AD')) vcf.B[[x]] <- as.character(vcf.B[[x]])
    
    vcf.B$GT[vcf.B$GT == '0' & vcf.A$CHROM == 'chrX'] <- '0/0'
    vcf.B$GT[vcf.B$GT == '1' & vcf.A$CHROM == 'chrX'] <- '1/1'
    vcf.B$GT <- gsub('|', '/', vcf.B$GT, fixed = T)
    vcf.B$DP[is.na(vcf.B$DP)] <- 0
    vcf.B$AF <- round(as.numeric(sapply(strsplit(vcf.B$AD, ','), function(x) x[2])) / vcf.B$DP, 3)
    vcf.B$AD <- NULL
    
    vcf <- cbind(vcf.A[which(hard.mask),], vcf.B[which(hard.mask),])
    if (once){
      pos.out <- scales::comma(vcf$POS, accuracy = 1)
      once = F
    }
    vcf$POS.out <- pos.out
    
    soft.mask <- vcf$DP < args$min.seq.depth | !(vcf$GT %in% c('.', '0/0', '0/1', '1/1'))
    vcf$GT[soft.mask] <- '.'
    vcf$DP[soft.mask] <- 0
    vcf$AF[soft.mask] <- NaN

    to.trio.format <- function(x){
      homoref = x == '0/0'
      hetero = x == '0/1'
      homoalt = x == '1/1'
      na = !(x %in% c('0/0', '0/1', '1/1'))
      x[homoref] <- 0
      x[hetero] <- 1
      x[homoalt] <- 2
      x[na] <- NA
      return(as.numeric(x))
    }
    vcf$GT.trio <- to.trio.format(vcf$GT)
    vcfs[[sample]] <- vcf
  }
  return(vcfs)
}

vcfs <- load.samples(args)

###
# bins
###

bins.annot <- GRanges(seqnames = chr.order, ranges = IRanges(start = 1, width = unlist(chr.lengths)))
bins.annot <- GRanges(as.data.frame(slidingWindows(bins.annot, width = args$bin.size * 1e6, step = args$bin.size * 1e6)))
bins.sample <- as.data.frame(bins.annot)[,1:3] ; colnames(bins.sample) <- c('CHROM', 'START', 'END')
bins.sample$MIDS <- round(bins.sample$START + (bins.sample$END - bins.sample$START + 1)/2)
bins <- list()
for (sample in c(args$singles, args$trio)){
  bins[[sample]] <- bins.sample
  bins[[sample]]$HOMO <- NA
  bins[[sample]]$HETERO <- NA
  
  gr <- GRanges(seqnames = vcfs[[sample]]$CHROM, ranges = IRanges(start = vcfs[[sample]]$POS, width = 1))
  hits <- nearest(bins.annot, gr, select = 'all')
  hits <- split(hits@to, hits@from)
  
  bins[[sample]]$HOMO <- sapply(hits, function(x) length(which(vcfs[[sample]]$GT[x] %in% c('1/1'))))
  bins[[sample]]$HETERO <- sapply(hits, function(x) length(which(vcfs[[sample]]$GT[x] %in% c('0/1'))))
  bins[[sample]]$DIFF <- bins[[sample]]$HOMO - bins[[sample]]$HETERO
  bins[[sample]]$COUNT <- bins[[sample]]$HOMO + bins[[sample]]$HETERO
  
  tot.homo <- sum(bins[[sample]]$HOMO)
  tot.hetero <- sum(bins[[sample]]$HETERO)
  
  fishers <- rep(NA, nrow(bins[[sample]]))
  for (i in 1:nrow(bins[[sample]])){
    a <- bins[[sample]]$HETERO[i]
    b <- bins[[sample]]$HOMO[i]
    mask <- bins[[sample]]$CHROM != bins[[sample]]$CHROM[i] & bins[[sample]]$CHROM %in% chr.order[1:22]
    c <- sum(bins[[sample]]$HETERO[mask]) - a
    d <- sum(bins[[sample]]$HOMO[mask]) - b
    if (a + b > args$min.support){
      fishers[i] <- fisher.test(matrix(c(a,b,c,d), ncol = 2), alternative = 'less')$p.value
    }
  }
  bins[[sample]]$FISHER <- round(fishers, 4)
  bins[[sample]]$DEVIANT <- bins[[sample]]$FISHER < .05
  bins[[sample]]$LOH <- NA
  
  step = floor(args$successive.bins.required/2)
  
  for (i in 1:nrow(bins[[sample]])){
    if (is.na(bins[[sample]]$DEVIANT[i])) next
    chr <- as.character(bins[[sample]]$CHROM[i])
    s <- max(which(bins[[sample]]$CHROM == chr)[1],
             i - floor(args$successive.bins.required / 2))
    e <- min(which(bins[[sample]]$CHROM == chr)[length(which(bins[[sample]]$CHROM == chr))],
             i + floor(args$successive.bins.required / 2))
    w <- s:e
    if (mean(is.na(bins[[sample]]$DEVIANT[w])) < .5){ # more than 50% NA -> NA
      if (all(bins[[sample]]$DEVIANT[w], na.rm = T)){
        bins[[sample]]$LOH[i] <- T
      } else {
        bins[[sample]]$LOH[i] <- F
      }
    }
  }
}

###
# bins.UDP.trio
###

if (length(args$trio)){
  geno.trio <- matrix(nrow = 3, ncol = nrow(vcfs[[1]]) + 2)
  for (i in 1:nrow(vcfs[[1]])){
    out.pos <- i + 2
    for (j in 1:3){
      geno.trio[j,out.pos] <- vcfs[[args$trio[j]]]$GT.trio[i]
    }
  }
  
  geno.trio[,1] <- 1; geno.trio[,2] <- 1:3
  geno.trio <- as.data.frame(geno.trio)
  colnames(geno.trio)[1:2] <- c('famid', 'pid')
  geno.trio$famid <- as.character(geno.trio$famid)
  geno.trio$pid <- as.character(geno.trio$pid)
  
  f = file()
  sink(file=f) ## silence output
  errors.i <- invisible(trio.check(geno.trio, is.linkage = F, replace = F)$errors$snp)
  sink() ## undo silencing
  close(f)
  
  GT.error <- t(geno.trio[,errors.i+2])
  GT.error.ref <- GT.error[which(GT.error[,3] == 0 & GT.error[,1] != GT.error[,2]),]
  GT.error.alt <- GT.error[which(GT.error[,3] == 2 & GT.error[,1] != GT.error[,2]),]
  
  father.support <- rownames(GT.error.ref)[which(GT.error.ref[,1] %in% c(0, 1))]
  mother.support <- rownames(GT.error.ref)[which(GT.error.ref[,2] %in% c(0, 1))]
  father.support <- c(father.support, rownames(GT.error.alt)[which(GT.error.alt[,1] %in% c(1, 2))])
  mother.support <- c(mother.support, rownames(GT.error.alt)[which(GT.error.alt[,2] %in% c(1, 2))])
  
  father.support.i <- sort(as.numeric(substring(father.support, 2)) - 2)
  mother.support.i <- sort(as.numeric(substring(mother.support, 2)) - 2)
  
  vcfs[[args$trio[3]]]$MEN.ER <- F
  vcfs[[args$trio[3]]]$MEN.ER[errors.i] <- T
  vcfs[[args$trio[3]]]$MEN.ER[is.na(vcfs[[args$trio[3]]]$GT.trio)] <- NA
  
  vcfs[[args$trio[3]]]$PAR <- NA
  vcfs[[args$trio[3]]]$PAR[father.support.i] <- 'fat'
  vcfs[[args$trio[3]]]$PAR[mother.support.i] <- 'mot'
  
  bins.trio.proband <- bins.sample
  gr <- GRanges(seqnames = vcfs[[args$trio[3]]]$CHROM, ranges = IRanges(start = vcfs[[args$trio[3]]]$POS, width = 1))
  hits <- nearest(bins.annot, gr, select = 'all')
  hits <- split(hits@to, hits@from)
  bins.trio.proband$MEN.ER <- sapply(hits, function(x) sum(vcfs[[args$trio[3]]]$MEN.ER[x], na.rm = T))
  bins.trio.proband$NO.MEN.ER <- sapply(hits, function(x) sum(!vcfs[[args$trio[3]]]$MEN.ER[x], na.rm = T))
  bins.trio.proband$FAT <- sapply(hits, function(x) sum(vcfs[[args$trio[3]]]$PAR[x] == 'fat', na.rm = T))
  bins.trio.proband$MOT <- sapply(hits, function(x) sum(vcfs[[args$trio[3]]]$PAR[x] == 'mot', na.rm = T))
  
  fishers <- rep(NA, nrow(bins.trio.proband))
  for (i in 1:nrow(bins.trio.proband)){
    a <- bins.trio.proband$NO.MEN.ER[i]
    b <- bins.trio.proband$MEN.ER[i]
    mask <- bins[[sample]]$CHROM != bins[[sample]]$CHROM[i] & bins[[sample]]$CHROM %in% chr.order[1:22]
    c <- sum(bins.trio.proband$NO.MEN.ER[mask]) - a
    d <- sum(bins.trio.proband$MEN.ER[mask]) - b
    if (a + b > args$min.support){
      fishers[i] <- fisher.test(matrix(c(a,b,c,d), ncol = 2), alternative = 'less')$p.value
    }
  }
  bins.trio.proband$FISHER <- round(fishers, 4)
  bins.trio.proband$DEVIANT <- bins.trio.proband$FISHER < .1
  bins.trio.proband$UDP <- NA
  
  step = floor(args$successive.bins.required/2)
  
  for (i in 1:nrow(bins.trio.proband)){
    if (is.na(bins.trio.proband$DEVIANT[i])) next
    chr <- as.character(bins.trio.proband$CHROM[i])
    s <- max(which(bins.trio.proband$CHROM == chr)[1],
             i - floor(args$successive.bins.required / 2))
    e <- min(which(bins.trio.proband$CHROM == chr)[length(which(bins.trio.proband$CHROM == chr))],
             i + floor(args$successive.bins.required / 2))
    w <- s:e
    if (mean(is.na(bins.trio.proband$DEVIANT[w])) < .5){ # more than 50% NA -> NA
      if (all(bins.trio.proband$DEVIANT[w], na.rm = T)){
        bins.trio.proband$UDP[i] <- T
      } else {
        bins.trio.proband$UDP[i] <- F
      }
    }
  }
}

###
# plots
###

zip <- function(...) unlist(mapply(list, ..., SIMPLIFY = FALSE))

plot.BAF <- function(sample, chr, trio = F, annot = 'this.is.annot'){ # uses l.matrix (see below)
  
  par(mar = c(.5,5,.5,0), mgp=c(0,.3,.3))
  plot.new()
  w <- which(vcfs[[sample]]$CHROM == chr)
  xlim = c(1, chr.lengths[[chr]])
  ylim = c(0, 100)
  plot(vcfs[[sample]]$POS[w], vcfs[[sample]]$AF[w] * 100, pch = 16, cex = .5, axes = F,
       ylab = '', xlab = '', ylim = ylim, xlim = xlim)
  rect(xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col = rgb(0,0,0,.05))
  axis(2, las=1, tcl=0.5, at = seq(ylim[1], ylim[2], length.out = 3))
  
  mtext('AF (%)', 4, -.7, cex=.8)
  
  text(0, 120, annot, cex = 1.5, adj = 0)
  
  par(mar = c(.5,1,.5,2))
  xhist <- hist(vcfs[[sample]]$AF[w] * 100, breaks = seq(0, 100, length.out = args$hist.bins), plot = F)
  xlim = c(0, max(xhist$counts))
  barplot(xhist$counts, axes = F, space = 0, horiz = T, xlab= "", ylab="", col = 'black')
  rect(xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col = rgb(0,0,0,.05))
  par(mgp=c(0,.3,.8))
  axis(3, las=1, tcl=0.5, at = xlim)
  par(mgp=c(0,.3,.3))
  mtext('count', 3, 1, cex=.8)
  
  par(mar = c(.5,5,.5,15))
  w <- which(bins[[sample]]$CHROM == chr)
  xlim = c(c(1, chr.lengths[[chr]]))
  ylim = c(0, max(bins[[sample]]$COUNT[w]) + (max(bins[[sample]]$COUNT[w]) %% 4))
  plot(bins[[sample]]$MIDS[w], bins[[sample]]$COUNT[w], type = 'l', axes = F,
       ylab = '', xlab = '', xlim = xlim, ylim = ylim, lwd = 0)
  rect(xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col = rgb(0,0,0,.05))
  polygon(zip(c(1, bins[[sample]]$START[w], chr.lengths[chr]),
              c(1, bins[[sample]]$END[w], chr.lengths[chr])),
          zip(c(0, bins[[sample]]$COUNT[w], 0),
              c(0, bins[[sample]]$COUNT[w], 0)), col = 'grey')
  axis(2, las=1, tcl=0.5, at = ylim)
  axis(4, las=1, tcl=0.5, at = ylim, labels = c('',''))
  text(chr.lengths[[chr]] + chr.lengths[[chr]] * .07, ylim[1] + diff(ylim) * .5, '#variants', cex=1.2, adj = 0)
  
  ylim = c(min(bins[[sample]]$DIFF[w]) + (min(bins[[sample]]$DIFF[w]) %% 4),
           max(bins[[sample]]$DIFF[w]) + (max(bins[[sample]]$DIFF[w]) %% 4))
  plot(zip(bins[[sample]]$START[w],
           bins[[sample]]$END[w]),
       zip(bins[[sample]]$DIFF[w],
           bins[[sample]]$DIFF[w]), type = 'l', axes = F,
       ylab = '', xlab = '', xlim = xlim, ylim = ylim)
  rect(xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col = rgb(0,0,0,.05))
  axis(2, las=1, tcl=0.5, at = ylim)
  axis(4, las=1, tcl=0.5, at = ylim, labels = c('', ''))
  text(chr.lengths[[chr]] + chr.lengths[[chr]] * .07, ylim[1] + diff(ylim) * .5, '#homo - #hetero', cex=1.2, adj = 0)
  
  ylim = c(0,1)
  plot(zip(bins[[sample]]$START[w],
           bins[[sample]]$END[w]),
       zip(bins[[sample]]$FISHER[w],
           bins[[sample]]$FISHER[w]), type = 'l', axes = F,
       ylab = '', xlab = '', xlim = xlim, lwd = 1, lty = 3, ylim = ylim)
  rect(xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col = rgb(0,0,0,.05))
  axis(2, las=1, tcl=0.5, at = seq(ylim[1], ylim[2], length.out = 3))
  axis(4, las=1, tcl=0.5, at = ylim, labels = c('', ''))
  text(chr.lengths[[chr]] + chr.lengths[[chr]] * .07, .5, 'P-value LOH (Fisher)', cex=1.2, adj = 0)
  
  segments(bins[[sample]]$START[w], ifelse(bins[[sample]]$LOH[w], -.2, NA),
           bins[[sample]]$END[w], ifelse(bins[[sample]]$LOH[w], -.2, NA), lwd = 3, col = color.A)
  segments(bins[[sample]]$START[w], ifelse(!bins[[sample]]$LOH[w], -.2, NA),
           bins[[sample]]$END[w], ifelse(!bins[[sample]]$LOH[w], -.2, NA), lwd = 3, col = color.C)
  text(chr.lengths[[chr]] + chr.lengths[[chr]] * .01, -.28, 'LOH:', cex=1, adj = 0, font = 2)
  text(chr.lengths[[chr]] + chr.lengths[[chr]] * .08, -.28, 'significant', cex=1, adj = 0, col = color.A)
  text(chr.lengths[[chr]] + chr.lengths[[chr]] * .2, -.28, 'not significant', cex=1, adj = 0, col = color.C)
  
  if (trio){
    
    ylim = c(0, max(c(bins.trio.proband$MEN.ER[w], 1)))
    plot(zip(bins.trio.proband$START[w],
             bins.trio.proband$END[w]),
         zip(bins.trio.proband$MEN.ER[w],
             bins.trio.proband$MEN.ER[w]), type = 'l', axes = F,
         ylab = '', xlab = '', xlim = xlim, lwd = 1, ylim = ylim)
    rect(xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col = rgb(0,0,0,.05))
    axis(2, las=1, tcl=0.5, at = ylim)
    axis(4, las=1, tcl=0.5, at = ylim, labels = c('', ''))
    text(chr.lengths[[chr]] + chr.lengths[[chr]] * .07, ylim[1] + diff(ylim) * .5, '#MEs', cex=1.2, adj = 0)
    
    ylim = c(0,1)
    plot(zip(bins.trio.proband$START[w],
             bins.trio.proband$END[w]),
         zip(bins.trio.proband$FISHER[w],
             bins.trio.proband$FISHER[w]), type = 'l', axes = F,
         ylab = '', xlab = '', xlim = xlim, lwd = 1, lty = 3, ylim = ylim)
    rect(xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col = rgb(0,0,0,.05))
    axis(2, las=1, tcl=0.5, at = seq(ylim[1], ylim[2], length.out = 3))
    axis(4, las=1, tcl=0.5, at = ylim, labels = c('', ''))
    text(chr.lengths[[chr]] + chr.lengths[[chr]] * .07, .5, 'P-value ME (Fisher)', cex=1.2, adj = 0)
    
    segments(bins.trio.proband$START[w], ifelse(bins.trio.proband$UDP[w], -.2, NA),
             bins.trio.proband$END[w], ifelse(bins.trio.proband$UDP[w], -.2, NA), lwd = 3, col = color.A)
    segments(bins.trio.proband$START[w], ifelse(!bins.trio.proband$UDP[w], -.2, NA),
             bins.trio.proband$END[w], ifelse(!bins.trio.proband$UDP[w], -.2, NA), lwd = 3, col = color.C)
    text(chr.lengths[[chr]] + chr.lengths[[chr]] * .01, -.28, 'ME:', cex=1, adj = 0, font = 2)
    text(chr.lengths[[chr]] + chr.lengths[[chr]] * .08, -.28, 'significant', cex=1, adj = 0, col = color.A)
    text(chr.lengths[[chr]] + chr.lengths[[chr]] * .2, -.28, 'not significant', cex=1, adj = 0, col = color.C)
    
    ylim = c(0, max(c(max(bins.trio.proband$FAT[w]), max(bins.trio.proband$MOT[w]), 1)))
    plot(zip(bins.trio.proband$START[w],
             bins.trio.proband$END[w]),
         zip(bins.trio.proband$FAT[w],
             bins.trio.proband$FAT[w]), type = 'l', axes = F,
         ylab = '', xlab = '', xlim = xlim, lwd = 1, lty = 1, ylim = ylim, col = color.A)
    par(new = T)
    plot(zip(bins.trio.proband$START[w],
             bins.trio.proband$END[w]),
         zip(bins.trio.proband$MOT[w],
             bins.trio.proband$MOT[w]), type = 'l', axes = F,
         ylab = '', xlab = '', xlim = xlim, lwd = 1, lty = 1, ylim = ylim, col = color.C)
    rect(xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col = rgb(0,0,0,.05))
    
    axis(2, las=1, tcl=0.5, at = ylim)
    axis(4, las=1, tcl=0.5, at = ylim, labels = c('', ''))
    text(chr.lengths[[chr]] + chr.lengths[[chr]] * .07, ylim[1] + diff(ylim) * .5, '#parent supporting MEs', cex=1.2, adj = 0)
    text(chr.lengths[[chr]] + chr.lengths[[chr]] * .14, ylim[1] + diff(ylim) * .1, 'father', cex=1, adj = 0, col = color.A)
    text(chr.lengths[[chr]] + chr.lengths[[chr]] * .22, ylim[1] + diff(ylim) * .1, 'mother', cex=1, adj = 0, col = color.C)
  }
  par(mgp=c(0,1,1.5))
  axis(1, las=1, tcl=0.5)
  plot.new()
}

get.sub.matrix <- function(n.sample, nrow = 7){
  l.matrix <- matrix(nrow = nrow, ncol = 4)
  l.matrix[1,] <- rep(1,4) + nrow * (n.sample - 1)
  l.matrix[2,] <- c(2,2,2,3) + nrow * (n.sample - 1)
  l.matrix[3,] <- c(2,2,2,3) + nrow * (n.sample - 1)
  l.matrix[4,] <- rep(4,4) + nrow * (n.sample - 1)
  l.matrix[5,] <- rep(5,4) + nrow * (n.sample - 1)
  l.matrix[6,] <- rep(6,4) + nrow * (n.sample - 1)
  l.matrix[7,] <- rep(7,4) + nrow * (n.sample - 1)
  return(l.matrix)
}
get.l.matrix <- function(n.samples){
  l.matrix <- matrix(nrow = 0, ncol = 4)
  for (n in 1:n.samples){
    l.matrix <- rbind(l.matrix, get.sub.matrix(n))
  }
  l.matrix <- rbind(l.matrix, rep(22, 4))
  l.matrix <- rbind(l.matrix, rep(23, 4))
  l.matrix <- rbind(l.matrix, rep(24, 4))
  return(l.matrix)
}

if (length(args$trio)){
  dir.create(paste0(args$output.dir, '/trio'), showWarnings = F, recursive = T)
  for (chr in chr.order){
    png(paste0(args$output.dir, '/trio/', chr, '.png'), width=8,height=10,units='in',res=512/2)
    par(oma = c(2,0,0,1), xpd = NA)
    layout(get.l.matrix(3))
    plot.BAF(args$trio[1], chr, annot = paste0('father: ', args$trio[1]))
    plot.BAF(args$trio[2], chr, annot = paste0('mother: ', args$trio[2]))
    plot.BAF(args$trio[3], chr, trio = T, annot = paste0('proband: ', args$trio[3]))
    dev.off()
  }
} else {
  for (sample in args$singles){
    dir.create(paste0(args$output.dir, '/', sample), showWarnings = F, recursive = T)
    for (group in split(chr.order, ceiling(seq_along(chr.order)/3))){
      png(paste0(args$output.dir, '/', sample, '/', paste0(group, collapse = '-'), '.png'),
          width=8,height=10,units='in',res=512/2)
      par(oma = c(2,0,0,0), xpd = NA)
      layout(get.l.matrix(3))
      for (i in 1:length(group)){
        plot.BAF(sample, group[i], annot = paste0(sample, ': ', group[i]))
      }
      dev.off()
    }
  }
}

###
# tables
###

for (sample in c(args$singles, args$trio)){
  full.output <- vcfs[[sample]][vcfs[[sample]]$DP != 0,][,c(1:2,4:8)]
  bin.output <- bins[[sample]][,-10]
  colnames(bin.output)[9] <- c('LOH.FISHER')
  if (length(args$trio)){
    if (sample == args$trio[3]){
      full.output <- vcfs[[sample]][vcfs[[sample]]$DP != 0,][,c(1:2, 4:8, 11, 12)]
      bin.output <- cbind(bin.output, bins.trio.proband[,c(5:9, 11)])
      colnames(bin.output)[15] <- c('UDP.FISHER')
    }
  }
  write.table(full.output, paste0(args$output.dir, '/', sample, '_pos.txt'), sep = '\t', quote = F, row.names = F, col.names = T)
  write.table(bin.output, paste0(args$output.dir, '/', sample, '_bin.txt'), sep = '\t', quote = F, row.names = F, col.names = T)
}



