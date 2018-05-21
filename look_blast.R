suppressMessages(library("data.table"))
suppressMessages(library(GenomicRanges))

max_dist <- 20000

f <- commandArgs(trailingOnly=TRUE)[1]
t <- fread(f, head = F)
colnames(t) <- c("qid", "sid", "identity", "length", "mismatches", "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
q <- strsplit(t$qid, "\\|")
q <- (do.call(rbind, q))
tab <- t
tab[,c("chr", "start", "end"):=.(q[,4], as.numeric(q[,5]),as.numeric(q[,6]))]
tab.orig <- tab


t <- (t[,.SD[which(bitscore==max(bitscore))], by=c("sid","qid")])
t <- (t[,.SD[which(bitscore==max(bitscore))], by="qid"])
ami.genes.count <- length(unique(t$sid))

fname <- gsub(".out", "", f)
fname <- gsub("blast/", "", fname)
tab <- t
q <- strsplit(t$qid, "\\|")
q <- (do.call(rbind, q))
tab[,c("chr", "start", "end"):=.(q[,4], as.numeric(q[,5]),as.numeric(q[,6]))]
tab <- tab[order(tab$chr,tab$start),]

gr <- makeGRangesFromDataFrame(tab[,c("chr", "start", "end","sid"), with=FALSE], keep.extra.columns=TRUE)
ok.dist.names <- list()
q <- 1
gr.ok <- GRanges()

good.start = -1
good.length = 1
max.good.length = 0
max.start = -1
max.end = -1

for(i in 1:(length(gr)-1)){
  .dist <- distance(gr[i],gr[i+1])
  if(is.na(.dist) | .dist > max_dist){
    if(good.length > max.good.length){
      max.good.length = good.length
      max.start = good.start
      max.end = i
    }
    good.start = -1 
    good.length = 1
    next
  }
  
  if(.dist < max_dist & good.start < 0) { 
    good.start = i
  }
  if(.dist < max_dist) { good.length = good.length + 1 }
}

cluster.range <- GRanges(seqnames = seqnames(gr[max.start]), IRanges(start(gr[max.start]), end(gr[max.end])))
tab.in.cluster <- tab[as.data.frame(findOverlaps(cluster.range, gr))[,2],]
write.table(tab.in.cluster, file=paste0("/data7a/bio/Shariki/PAPER2/detectAmiManyBac/blast_table/",fname,".blast_res"), quote=F, row.names=F, sep="\t", col.names = F)

names. <- substr(tab.in.cluster$sid,1,4)
essential.genes.present <- sum(c("AmiA","AmiB","AmiI","AmiK","AmiL", "AmiM") %in% names.)
all.genes <- paste0("Ami", LETTERS[1:15])
genes.present <- sum(all.genes %in% names.)

gr.orig <- makeGRangesFromDataFrame(tab.orig[,c("chr", "start", "end","sid"), with=FALSE], keep.extra.columns=TRUE)
tab.orig.in.cluster <- tab.orig[as.data.frame(findOverlaps(cluster.range, gr.orig))[,2],]
to <- tab.orig[,c("sid", "sstart", "send"), with=FALSE]
colnames(to) <- c("seqname","start","end")
cluster.ranges <- makeGRangesFromDataFrame(to)
cluster.ranges <- reduce(cluster.ranges)
cluster.ranges <- data.frame(cluster.ranges)
cluster.ranges$seqnames <- as.character(cluster.ranges$seqnames)
cluster.ranges.agg <- aggregate(width ~ seqnames, cluster.ranges,sum)
write.table(cluster.ranges.agg, file=paste0("/data7a/bio/Shariki/PAPER2/detectAmiManyBac/cluster_gene_blastlength/",fname,".bl"), quote=F, row.names=F, sep="\t", col.names = F)


cat(c(fname, genes.present, essential.genes.present, sum_bitscore=sum(tab.in.cluster$bitscore),sum_length=sum(cluster.ranges.agg$width)), file=paste0("/data1/bio/Shariki/allbac/blast_stat/",fname,".blast_stat4"), sep = '\t')

write.table(c(fname, names.), file=paste0("/data1/bio/Shariki/allbac/blast_stat/names/",fname,".names"), quote=F, row.names=F, sep="\t", col.names = F)
