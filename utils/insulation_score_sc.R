# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("viridis")
  library("ggpointdensity")
  library("GenomicRanges")
  library("ggpubr")
})

# result folder
result_folder = "../results/hi-c/"

# Hi-C boundaries
hic_boundary_beds = "../data/Hi-C/bed/"
list.files(input, full.names = TRUE)

# scatterplot function
# it visualize the log2 insulation scores of the determined boundaries
# you feed two bed files (full path) and two axis titles
boundary_scatter = function(boundaries_1,
                            boundaries_2,
                            x_title, 
                            y_title
                            ) {
  
  boundaries_1 = fread(boundaries_1)
  boundaries_2 = fread(boundaries_2)
  
  boundaries_1 = GRanges(
    seqnames = boundaries_1$seqnames,
    ranges = IRanges(
      start = boundaries_1$start,
      end = boundaries_1$end,
      insulation_score = boundaries_1$log2_insulation_score_200000
    )
  )
  
  boundaries_2 = GRanges(
    seqnames = boundaries_2$seqnames,
    ranges = IRanges(
      start = boundaries_2$start,
      end = boundaries_2$end,
      insulation_score = boundaries_2$log2_insulation_score_200000
    )
  )
  
  ol = findOverlaps(boundaries_1,
                    boundaries_2,
                    type = "any",
                    ignore.strand = FALSE)
  
  if(length(ol) == 0) {
    return(NULL)
  }
  
  
  boundaries_1 = boundaries_1[queryHits(ol)]
  boundaries_2 = boundaries_2[subjectHits(ol)]
  
  boundaries_1 = as.data.frame(boundaries_1)
  boundaries_2 = as.data.frame(boundaries_2)
  
  vis = tibble(boundaries_1 = boundaries_1$insulation_score, boundaries_2 = boundaries_2$insulation_score)
  
  
  sc = ggplot(data = vis, mapping = aes(x = boundaries_1, y = boundaries_2)) +
    geom_pointdensity(show.legend = FALSE) +
    scale_colour_gradient(low = "#f03b20",
                          high = "yellow",
                          breaks = c(500, 1000, 1500, 2000)) +
    labs(colour = " ") +
    ylab(y_title) +
    xlab(x_title) +
    theme_minimal() +
    xlim(-10, 3) +
    ylim(-10, 3) +
    theme(
      axis.line = element_line(size = 0.5, colour = "black"),
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(
        size = 13,
        color = "black",
        angle = 0,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(size = 13, color = "black")
    ) +
    stat_cor(method = "pearson", label.x = -9.0, label.y = 3, size = 6)
  
  return(print(sc))
  
  
}

boundary_scatter(boundaries_1 = "../data/Hi-C/bed/common_fastq_merge_200kb.bed", 
                 boundaries_2 = "../data/Hi-C/bed/NT_fastq_merge_bd_200kb.bed",
                 y_title = "common", x_title = "non-treated merge")

boundaries = list.files(hic_boundary_beds, full.names = TRUE, pattern = "*_200kb.bed")
labels = c("common", "non-treated merge", "non-treated only", "TMPyP4", "TMPyP4 only")
sc_list = list()
counter = 0
for (i in 1:length(boundaries)) {
  for (j in 1:length(boundaries)) {
    counter = ifelse(i %in% 1:length(boundaries), counter + 1, counter)
    plot = boundary_scatter(boundaries_1 = boundaries[i],
                            boundaries_2 = boundaries[j],
                            x_title = labels[i],
                            y_title = labels[j])
    
    sc_list[[as.character(counter)]] = plot
  }
}
sc_list[sapply(sc_list, is.null)] = NULL
ggarrange(plotlist = sc_list)

ggsave(
  plot = last_plot(),
  glue("{result_folder}boundary_scatters.pdf"),
  width = 14,
  height = 12
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}boundary_scatters.png"),
  width = 14,
  height = 12,
  dpi = 300
)
  