library("pqsfinder")
library("data.table")
library("rtracklayer")

# PQS finder algorithm
# source: http://127.0.0.1:10220/library/pqsfinder/doc/pqsfinder.html#g-quadruplex-detection

# needs a fasta: https://hgdownload.soe.ucsc.edu/downloads.html
mm10 = readDNAStringSet(file = "../data/ucsc_mm10.fa", format = "fasta")

# pqsfinder's function
# IT TAKES MORE THAN 1 HOUR!!!
pqs = lapply(mm10, pqsfinder)

# create data frame
pqs_gr = lapply(pqs, as, Class = "GRanges")
pqs_gr_unlisted = unlist(pqs_gr)
pqs_gr_dfs = lapply(pqs_gr_unlisted, as.data.frame)
pqs_gr_df = rbindlist(pqs_gr_dfs)
pqs_gr_df = pqs_gr_df[,c("seqnames", "start", "end", "width", "strand", "score")]

# export
save(pqs_gr_df, file = "../data/mm10_PQS.Rda")
write_tsv(pqs_gr_df, "../data/mm10_PQS.tsv")

# export in bigwig format
pqs_gr_df_gr = makeGRangesFromDataFrame(pqs_gr_df,
                         keep.extra.columns = TRUE,
                         ignore.strand = FALSE,
                         seqinfo = NULL,
                         seqnames.field = c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field = "start",
                         end.field = c("end", "stop"),
                         strand.field = "strand",
                         starts.in.df.are.0based = FALSE)
export.bw(coverage(pqs_gr_df_gr), "../data/mm10_PQS.bw")
