library("rtracklayer")
library("GenomicRanges")

# Rui's script to convert bed files to bigwigs (e.g. PQS analysis)
# readPQS function: specific to Rui's G4 site predicting script
readPQS <- function(file_path)
{
  tmp <- read.table(file_path, 
                    header = F, sep = '\t',
                    col.names = c('seqnames', 'start', 'end', 'name', 'length', 'strand', 'seq')
  )
  
  tmp$seqnames <- gsub('(.*)\\sN.*', '\\1', tmp$seqnames)
  tmp <- tmp[tmp$seqnames %in% paste0('chr', c(1:100, 'X', 'Y')), ]
  makeGRangesFromDataFrame(tmp)
}

# converting bed to bigwig (mm10 or hg19)
bed_to_bw <- function(file_path, .name, genome) {
  if (length(file_path) > 1) {
    gr <- Reduce("c", sapply(file_path, readPQS))
  } else {
    gr <- readPQS(file_path)
  }
  
  export.bw(coverage(gr), paste0(.name, "_stack.bw"))
  
  gr <- reduce(gr, ignore.strand = T)
  gr$score <- 1
  if (genome == "hg19") {
    seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[seqlevels(gr)]
  } 
  if (genome == "mm10") {
    seqlengths(gr) <- seqlengths(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)[seqlevels(gr)]
  }
  export.bw(gr, paste0(.name, "_binary.bw"))
}

bed_to_bw("../data/mm10_canPQS-regex.bed", "../data/mm10_canPQS-regex", genome = "mm10")
