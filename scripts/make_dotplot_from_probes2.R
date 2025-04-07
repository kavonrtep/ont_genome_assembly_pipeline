#!/usr/bin/env Rscript

library(optparse)
options_list <- list(
  make_option(c("-q","--query"), action = "store", type = "character",
              help = "query is blast output in tabular format", default = NULL),
  make_option(c("-s", "--subject"), action = "store", type = "character",
              help = "subject blast output in tabular format", default = NULL),
  make_option(c("-Q", "--query_chromosome_sizes"), action = "store", type = "character",
              help = "file with chromosome sizes for query(.fai format)", default = NULL),
  make_option(c("-S", "--subject_chromosome_sizes"), action = "store", type = "character",
              help = "file with chromosome sizes for subject(.fai format)", default = NULL),
  make_option(c("-o", "--output"), action = "store", type = "character",
              help = "output png file", default = NULL)
)

opt_parser <- OptionParser(option_list = options_list, description = "Creates dotplot from two blast outputs")
opt <- parse_args(opt_parser, args = commandArgs(TRUE))

# For testing
if (FALSE){
  opt <- list(
    query = "path/to/query.blast_out",
    subject = "path/to/subject.blast_out",
    query_chromosome_sizes = "path/to/query.fai",
    subject_chromosome_sizes = "path/to/subject.fai",
    output = "test.png"
  )
}

# Read BLAST tables (format 6)
query <- read.table(file = opt$query, sep="\t", header = FALSE,
                    col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                  "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
subject <- read.table(file = opt$subject, sep="\t", header = FALSE,
                      col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                    "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

# Read chromosome sizes in .fai format
chromosome_sizes_query <- read.table(file = opt$query_chromosome_sizes, sep="\t", header = FALSE,
                                     col.names = c('qseqid', 'qlen', 'V3','V4', 'V5'))
chromosome_sizes_subject <- read.table(file = opt$subject_chromosome_sizes, sep="\t", header = FALSE,
                                       col.names = c('sseqid', 'slen', 'V3','V4', 'V5'))

# --- NEW: Reorder subject contigs by length (largest first) ---
chromosome_sizes_subject <- chromosome_sizes_subject[order(as.numeric(chromosome_sizes_subject$slen), decreasing = TRUE), ]
# Calculate cumulative size for subject
chromosome_sizes_subject$cumulative_slen <- cumsum(as.numeric(chromosome_sizes_subject$slen))
chromosome_sizes_subject$cumulative_start <- c(1, head(chromosome_sizes_subject$cumulative_slen, -1))

# Determine strands from BLAST hit coordinates
query$qstrand <- ifelse(query$sstart < query$send, "plus", "minus")
subject$sstrand <- ifelse(subject$sstart < subject$send, "plus", "minus")

# Merge query and subject tables by qseqid (probe name)
QS <- merge(query, subject, by = "qseqid", suffixes = c("_query", "_subject"))

# (Optionally inspect strand table)
table(QS$qstrand, QS$sstrand)

# For the merged table, we keep the subject strand as the representative strand
QS$strand <- ifelse(QS$qstrand == QS$sstrand, "plus", "minus")

new_colnames <- c(
  sseqid_query = "qseqid",
  sseqid_subject = "sseqid",
  sstart_query = "qstart",
  send_query = "qend",
  sstart_subject = "sstart",
  send_subject = "send",
  strand = "sstrand"
)

blast_df <- QS[, names(new_colnames)]
colnames(blast_df) <- new_colnames

# --- Compute subject relative coordinates using the reordered subject ---
blast_df$sstart_relative <- blast_df$sstart + chromosome_sizes_subject$cumulative_start[match(blast_df$sseqid, chromosome_sizes_subject$sseqid)]
blast_df$send_relative <- blast_df$send + chromosome_sizes_subject$cumulative_start[match(blast_df$sseqid, chromosome_sizes_subject$sseqid)]

# --- NEW: Reorder query contigs based on median subject coordinate ---
# Compute a midpoint for each BLAST hit in subject coordinates

blast_df$subject_mid <- (blast_df$sstart_relative + blast_df$send_relative) / 2
# Calculate the median subject coordinate for each query contig
median_subject <- aggregate(subject_mid ~ qseqid, data = blast_df, FUN = median)
# Merge the median info into the query chromosome sizes table
chromosome_sizes_query <- merge(chromosome_sizes_query, median_subject, by = "qseqid", all.x = TRUE)
# Order query contigs by the median subject coordinate (ascending)
chromosome_sizes_query <- chromosome_sizes_query[order(chromosome_sizes_query$subject_mid), ]
# Recalculate cumulative sizes for query based on new order
chromosome_sizes_query$cumulative_qlen <- cumsum(as.numeric(chromosome_sizes_query$qlen))
chromosome_sizes_query$cumulative_start <- c(1, head(chromosome_sizes_query$cumulative_qlen, -1))

# --- Update factor levels for plotting ---
blast_df$qseqid <- factor(blast_df$qseqid, levels = chromosome_sizes_query$qseqid)
blast_df$sseqid <- factor(blast_df$sseqid, levels = chromosome_sizes_subject$sseqid)

# Now compute query relative coordinates using the new query ordering
blast_df$qstart_relative <- blast_df$qstart + chromosome_sizes_query$cumulative_start[match(blast_df$qseqid, chromosome_sizes_query$qseqid)]
blast_df$qend_relative <- blast_df$qend + chromosome_sizes_query$cumulative_start[match(blast_df$qseqid, chromosome_sizes_query$qseqid)]

# --- Plotting ---
# Calculate maximum coordinates for plot limits
xmax <- max(chromosome_sizes_query$cumulative_qlen)
ymax <- max(chromosome_sizes_subject$cumulative_slen)

# Define plot dimensions (height fixed; width proportional to subject genome size)
H <- 10
M <- 1
W <- round(H * xmax / ymax)
query_label <- basename(opt$query)
subject_label <- basename(opt$subject)

warning("Creating plot")
warning("png dimension 1: ", W + M * 2)
warning("png dimension 2: ", H + M * 2)
width <- W + M * 2
height <- H + M * 2
if (width < 5 | height < 5){
  width <- width * 2
  height <- height * 2
}

tryCatch({
  png(opt$output, width = width, height = height, res = 600, units = "in")
}, error = function(e){
  message("Error creating png file")
  message(e)
  png(opt$output, width = 12, height = 122, res = 600, units = "in")
})

par(mai = rep(M, 4))
plot(1, 1, xlim = c(1, xmax), ylim = c(1, ymax),
     xlab = query_label, ylab = subject_label, type = "n", axes = FALSE)

# Draw boundaries for each chromosome (query and subject)
par(lwd = 1)
abline(v = c(1, chromosome_sizes_query$cumulative_qlen), lty = 2, col = "grey")
abline(h = c(1, chromosome_sizes_subject$cumulative_slen), lty = 2, col = "grey")

# Add chromosome labels
axis(side = 1, at = chromosome_sizes_query$cumulative_qlen - chromosome_sizes_query$qlen / 2,
     labels = chromosome_sizes_query$qseqid, tick = FALSE, cex.axis = 0.5, las = 2)
axis(side = 2, at = chromosome_sizes_subject$cumulative_slen - chromosome_sizes_subject$slen / 2,
     labels = chromosome_sizes_subject$sseqid, tick = FALSE, cex.axis = 0.5, las = 2)

# Define color scheme for hits based on strand (blue for plus, red for minus)
color_scheme <- ifelse(blast_df$sstrand == "plus", "blue", "red")

# Draw segments for each BLAST hit (using relative coordinates)
segments(blast_df$qstart_relative, blast_df$sstart_relative,
         blast_df$qend_relative, blast_df$send_relative,
         col = color_scheme)

dev.off()




# --- Export data in PAF format ---
blast_df$qlen <- chromosome_sizes_query$qlen[match(blast_df$qseqid, chromosome_sizes_query$qseqid)]
blast_df$slen <- chromosome_sizes_subject$slen[match(blast_df$sseqid, chromosome_sizes_subject$sseqid)]
blast_df$length <- abs(blast_df$qend - blast_df$qstart)
blast_df$strand <- ifelse(blast_df$sstrand == "plus", "+", "-")
blast_df$mapq <- 255

blast_df$cumulative_qlen <- chromosome_sizes_query$cumulative_qlen[match(blast_df$qseqid, chromosome_sizes_query$qseqid)]
blast_df$cumulative_slen <- chromosome_sizes_subject$cumulative_slen[match(blast_df$sseqid, chromosome_sizes_subject$sseqid)]
blast_df$cumulative_start_qlen <- chromosome_sizes_query$cumulative_start[match(blast_df$qseqid, chromosome_sizes_query$qseqid)]
blast_df$cumulative_strat_slen <- chromosome_sizes_subject$cumulative_start[match(blast_df$sseqid, chromosome_sizes_subject$sseqid)]


if (is.null(opt$output_paf)) {
  opt$output_paf <- paste0(opt$output, ".paf")
}

PAF <- blast_df[, c("qseqid", "qlen", "qstart", "qend", "strand",
                      "sseqid", "slen", "sstart", "send", "length", "length", "mapq")]
colnames(PAF) <- c("qseqid", "qlen", "qstart", "qend", "strand",
                   "sseqid", "slen", "sstart", "send", "length", "m", "mapq")
write.table(PAF, file = opt$output_paf, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

tbl <- table(blast_df$qseqid, color_scheme)

# prevalent orientation
rc_orientation_id <- rownames(tbl)[tbl[,2]>tbl[,1]]


blast_df$qstart_relative_new <- blast_df$cumulative_qlen - blast_df$qstart_relative + blast_df$cumulative_start_qlen
blast_df$qend_relative_new   <- blast_df$cumulative_qlen - blast_df$qend_relative + blast_df$cumulative_start_qlen


# export second image
second_output <- paste0(tools::file_path_sans_ext(opt$output), "_reversed.png")

png(second_output, width = width, height = height, res = 600, units = "in")
par(mai = rep(M, 4))
plot(1, 1, xlim = c(1, xmax), ylim = c(1, ymax),
     xlab = query_label, ylab = subject_label, type = "n", axes = FALSE)

# Draw boundaries for each chromosome (query and subject)
par(lwd = 1)
abline(v = c(1, chromosome_sizes_query$cumulative_qlen), lty = 2, col = "grey")
abline(h = c(1, chromosome_sizes_subject$cumulative_slen), lty = 2, col = "grey")

# Add chromosome labels
axis(side = 1, at = chromosome_sizes_query$cumulative_qlen - chromosome_sizes_query$qlen / 2,
     labels = chromosome_sizes_query$qseqid, tick = FALSE, cex.axis = 0.5, las = 2)
axis(side = 2, at = chromosome_sizes_subject$cumulative_slen - chromosome_sizes_subject$slen / 2,
     labels = chromosome_sizes_subject$sseqid, tick = FALSE, cex.axis = 0.5, las = 2)

# Define color scheme for hits based on strand (blue for plus, red for minus)
color_scheme <- ifelse(blast_df$sstrand == "plus", "blue", "red")

# Draw segments for each BLAST hit (using relative coordinates)
x0 <- ifelse(blast_df$qseqid %in% rc_orientation_id, blast_df$qstart_relative_new, blast_df$qstart_relative)
x1 <- ifelse(blast_df$qseqid %in% rc_orientation_id, blast_df$qend_relative_new, blast_df$qend_relative)
y0 <- blast_df$sstart_relative
y1 <- blast_df$send_relative

segments(x0, y0, x1, y1,
         col = color_scheme)

dev.off()



