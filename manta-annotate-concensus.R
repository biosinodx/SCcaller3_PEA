args = commandArgs(trailingOnly=TRUE)
sn=args[1]

anchoredBases=150

library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)
#' Simple SV type classifier
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX",
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))
}

#
vcf <- readVcf(paste(sn, "/results/variants/diploidSV.vcf", sep=''), "hg38")
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
	row.names=c("SIMPLE_TYPE"),
	Number=c("1"),
	Type=c("String"),
	Description=c("Simple event type annotation based purely on breakend position and orientation."))), "DataFrame"))

# for humans, filter chromosome
c1=seqnames(vcf)=="chr1" | seqnames(vcf)=="chr2" | seqnames(vcf)=="chr3" | seqnames(vcf)=="chr4" | seqnames(vcf)=="chr5"
c2=seqnames(vcf)=="chr6" | seqnames(vcf)=="chr7" | seqnames(vcf)=="chr8" | seqnames(vcf)=="chr9" | seqnames(vcf)=="chr10"
c3=seqnames(vcf)=="chr11" | seqnames(vcf)=="chr12" | seqnames(vcf)=="chr13" | seqnames(vcf)=="chr14" | seqnames(vcf)=="chr15"
c4=seqnames(vcf)=="chr16" | seqnames(vcf)=="chr17" | seqnames(vcf)=="chr18" | seqnames(vcf)=="chr19" | seqnames(vcf)=="chr20"
c5=seqnames(vcf)=="chr21" | seqnames(vcf)=="chr22" | seqnames(vcf)=="chrX" | seqnames(vcf)=="chrY"
chr=c1 | c2 | c3 | c4 | c5
vcf_filt=vcf[which(chr), ]

# filter quaility
vcf_filt <- vcf_filt[which(rowRanges(vcf_filt)$FILTER=="PASS"), ]

# 
gr <- breakpointRanges(vcf_filt); # unpaired breakpoints are automatically removed
svtype <- simpleEventType(gr)
info(vcf_filt)$SIMPLE_TYPE <- NA_character_
info(vcf_filt[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf_filt[gr$sourceId])$SVLEN <- gr$svLen

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("BSgenome")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library("BSgenome")
library("Biostrings")
hg38 <- getBSgenome("hg38")
gr_p = partner(gr)

rsv=vector()
ra=vector()
rb=vector()
for(i in 1:length(gr)){
	if(i%%100==0)	print(i)
	rsv[i] <- StructuralVariantAnnotation:::extractBreakpointSequence(c(gr[i], gr_p[i]), hg38, anchoredBases=anchoredBases)[1]
	tmp <- StructuralVariantAnnotation:::extractReferenceSequence(c(gr[i], gr_p[i]), hg38, anchoredBases=anchoredBases)
	ra[i] = tmp[1]
	rb[i] = tmp[2]
}

info(header(vcf_filt)) = unique(as(rbind(as.data.frame(info(header(vcf_filt))), data.frame(
	row.names=c("REFSV"),
	Number=c("1"),
	Type=c("String"),
	Description=c("REFSV."))), "DataFrame"))
info(vcf_filt)$REFSV <- NA_character_
info(vcf_filt[gr$sourceId])$REFSV <- rsv

info(header(vcf_filt)) = unique(as(rbind(as.data.frame(info(header(vcf_filt))), data.frame(
	row.names=c("REFA"),
	Number=c("1"),
	Type=c("String"),
	Description=c("REFA."))), "DataFrame"))
info(vcf_filt)$REFA <- NA_character_
info(vcf_filt[gr$sourceId])$REFA <- ra

info(header(vcf_filt)) = unique(as(rbind(as.data.frame(info(header(vcf_filt))), data.frame(
	row.names=c("REFB"),
	Number=c("1"),
	Type=c("String"),
	Description=c("REFB."))), "DataFrame"))
info(vcf_filt)$REFB <- NA_character_
info(vcf_filt[gr$sourceId])$REFB <- rb

info(header(vcf_filt)) = unique(as(rbind(as.data.frame(info(header(vcf_filt))), data.frame(
	row.names=c("CHR2"),
	Number=c("1"),
	Type=c("String"),
	Description=c("CHR2."))), "DataFrame"))
info(vcf_filt)$CHR2 <- NA_character_
info(vcf_filt[gr$sourceId])$CHR2 <- seqnames(gr_p)

info(header(vcf_filt)) = unique(as(rbind(as.data.frame(info(header(vcf_filt))), data.frame(
	row.names=c("END2"),
	Number=c("1"),
	Type=c("Integer"),
	Description=c("END2."))), "DataFrame"))
info(vcf_filt)$END2 <- NaN
info(vcf_filt[gr$sourceId])$END2 <- start(gr_p)

writeVcf(vcf_filt, paste(sn, ".manta.concensus.vcf", sep=''))

