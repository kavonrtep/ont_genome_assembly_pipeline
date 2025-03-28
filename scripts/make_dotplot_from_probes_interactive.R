#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)

# Define command-line options
options_list <- list(
  make_option(c("-p", "--paf"), action = "store", type = "character",
              help = "Input PAF file", default = NULL),
  make_option(c("-Q", "--query_fai"), action = "store", type = "character",
              help = "Query FAI file", default = NULL),
  make_option(c("-S", "--subject_fai"), action = "store", type = "character",
              help = "Subject FAI file", default = NULL),
  make_option(c("-o", "--static_output"), action = "store", type = "character",
              help = "Output static PNG file", default = NULL),
  make_option(c("-i", "--interactive_output"), action = "store", type = "character",
              help = "Output interactive HTML file", default = NULL)
)

opt_parser <- OptionParser(option_list = options_list,
                           description = "Creates static and interactive dotplots from a PAF file and FAI files")
opt <- parse_args(opt_parser, args = commandArgs(TRUE))

#---------------------------
# Read input files
#---------------------------
# Read PAF file (no header; standard 12 columns)
paf <- read.table(file = opt$paf, sep = "\t", header = FALSE,
                  col.names = c("qseqid", "qlen", "qstart", "qend", "strand",
                                "sseqid", "slen", "sstart", "send", "length", "m", "mapq"))

# Read FAI files for query and subject.
# The FAI files are assumed to be in the standard format: sequence, length, and other columns.
query_fai <- read.table(file = opt$query_fai, sep = "\t", header = FALSE,
                        col.names = c("qseqid", "qlen", "V3", "V4", "V5"))
subject_fai <- read.table(file = opt$subject_fai, sep = "\t", header = FALSE,
                          col.names = c("sseqid", "slen", "V3", "V4", "V5"))

#---------------------------
# Compute cumulative coordinates
#---------------------------
# For query FAI:
query_fai$cumulative_qlen <- cumsum(as.numeric(query_fai$qlen))
query_fai$cumulative_start <- c(1, head(query_fai$cumulative_qlen, -1))
# For subject FAI:
subject_fai$cumulative_slen <- cumsum(as.numeric(subject_fai$slen))
subject_fai$cumulative_start <- c(1, head(subject_fai$cumulative_slen, -1))

#---------------------------
# Compute global (relative) coordinates for PAF segments
#---------------------------
paf$qstart_global <- paf$qstart + query_fai$cumulative_start[match(paf$qseqid, query_fai$qseqid)]
paf$qend_global   <- paf$qend   + query_fai$cumulative_start[match(paf$qseqid, query_fai$qseqid)]
paf$sstart_global <- paf$sstart + subject_fai$cumulative_start[match(paf$sseqid, subject_fai$sseqid)]
paf$send_global   <- paf$send   + subject_fai$cumulative_start[match(paf$sseqid, subject_fai$sseqid)]

# Order factor levels so that chromosomes appear in the correct order in the plot.
paf$qseqid <- factor(paf$qseqid, levels = query_fai$qseqid)
paf$sseqid <- factor(paf$sseqid, levels = subject_fai$sseqid)
save.image("debug.RData")
#---------------------------
# Create dotplot using ggplot2
#---------------------------
p <- ggplot(paf) +
  geom_segment(aes(x = qstart_global, xend = qend_global,
                   y = sstart_global, yend = send_global,
                   color = strand,
                   text = paste0("Query: ", qseqid, " (", qstart, "-", qend, ")",
                                 "<br>Subject: ", sseqid, " (", sstart, "-", send, ")")
                   ), lineend="square", size = 0.5
  ) +
  geom_vline(xintercept = c(0, query_fai$cumulative_qlen), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = c(0, subject_fai$cumulative_slen), linetype = "dashed", color = "grey") +
  scale_x_continuous(name = basename(opt$query_fai),
                     breaks = query_fai$cumulative_qlen - query_fai$qlen/2,
                     labels = query_fai$qseqid,
                     limits = c(0, max(query_fai$cumulative_qlen))) +
  scale_y_continuous(name = basename(opt$subject_fai),
                     breaks = subject_fai$cumulative_slen - subject_fai$slen/2,
                     labels = subject_fai$sseqid,
                     limits = c(0, max(subject_fai$cumulative_slen))) +
  scale_color_manual(values = c("+" = "blue", "-" = "red"), name = "Strand") +
  theme_minimal()

#---------------------------
# Save static dotplot as PNG
#---------------------------
# Set plot dimensions: fix height at 10 inches and scale width by the genome length ratio.
H <- 10
W <- H * max(query_fai$cumulative_qlen) / max(subject_fai$cumulative_slen)
ggsave(filename = opt$static_output, plot = p, width = W, height = H, dpi = 900)

#---------------------------
# Save interactive dotplot as HTML
#---------------------------
p_interactive <- ggplotly(p, tooltip = "text")
saveWidget(p_interactive, file = opt$interactive_output, selfcontained = TRUE)
