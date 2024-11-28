#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# print(args)
# Check provided arguments ----
if (length(args)!=10) {
  stop("Please specify:\n
       [required] 1/ A comma-separated list of sample count files\n
       [required] 2/ A comma-separated list of input count files\n
       [required] 3/ A comma-separated list of sample names (mathcing the list of count files)\n
       [required] 4/ A comma-separated list of input names (mathcing the list of count files)\n
       [required] 5/ Function to be applied to input columns for filtering\n
       [required] 6/ Function to be applied to sample columns for filtering\n
       [required] 7/ Function to be applied to all columns for filtering\n
       [required] 8/ Pseudocount to be added to input and sample columns (which will be further normalized for sequencing depth) \n
       [required] 9/ Raw counts output file\n
       [required] 10/ Filtered counts output file\n")
}

require(data.table)

# Test arguments ----
# sample_files <- c("db/counts/ORFeome/p53_nutlin_dim_A549_rep1_lib200_counts.txt", "db/counts/ORFeome//p53_nutlin_dim_A549_rep2_lib200_counts.txt")
# input_files <- c("db/counts/ORFeome/p53_input_NA_A549_rep1_lib200_counts.txt", "db/counts/ORFeome//p53_input_NA_A549_rep2_lib200_counts.txt")
# sample_names <- c("p53_nutlin_dim_A549_rep1", "p53_nutlin_dim_A549_rep2")
# input_names <- c("p53_input_NA_A549_rep1", "p53_input_NA_A549_rep2")
# sample.cutoff.FUN <- function (x)  sum(x) >= 0
# input.cutoff.FUN <- function (x)  sum(x) >= 0
# row.cutoff.FUN <- function (x)  sum(x) >= 1
# pseudocount <- 0
# raw_output <- "db/FC_tables/ORFeome//p53_nutlin_dim_A549/p53_nutlin_dim_A549_raw_counts.txt"
# filtered_output <- "db/FC_tables/ORFeome//p53_nutlin_dim_A549/p53_nutlin_dim_A549_filtered_counts.txt"

# Parse arguments ----
sample_files <- unlist(tstrsplit(args[1], ","))
input_files <- unlist(tstrsplit(args[2], ","))
sample_names <- unlist(tstrsplit(args[3], ","))
input_names <- unlist(tstrsplit(args[4], ","))
sample.cutoff.FUN <- eval(parse(text= args[5]))
input.cutoff.FUN <- eval(parse(text= args[6]))
row.cutoff.FUN <- eval(parse(text= args[7]))
pseudocount <- as.numeric(args[8])
raw_output <- args[9]
filtered_output <- args[10]

# Import counts ----
counts <- lapply(c(sample_files, input_files), fread)
names(counts) <- c(sample_names, input_names)
counts <- rbindlist(counts, idcol = "sampleID")
counts <- dcast(counts, ID~sampleID, value.var = "count")
setnames(counts, "ID", "sgRNA")
counts[, sgRNA:= gsub(",", "_", sgRNA)]
counts[, `gene name`:= tstrsplit(sgRNA, "__", keep= 1)]
setcolorder(counts,
            c("sgRNA", "gene name", sample_names, input_names))

# Save raw counts file ----
fwrite(counts,
       raw_output,
       col.names = TRUE,
       row.names = FALSE,
       sep = "\t",
       quote= F,
       na = NA)
total <- nrow(counts)

# Apply user-defined cutoffs ----
sample.sel <- apply(counts[, sample_names, with= F], 1, sample.cutoff.FUN)
input.sel <- apply(counts[, input_names, with= F], 1, input.cutoff.FUN)
row.sel <- apply(counts[, c(input_names, sample_names), with= F], 1, row.cutoff.FUN)
counts <- counts[sample.sel & input.sel & row.sel]
fil <- nrow(counts)
print(paste0(formatC(fil, big.mark = ","), " / ", formatC(total, big.mark = ","), " reads (", round(fil/total*100, 1), "%) passed sample/input cutoffs"))

# Apply user-defined pseudocount ----
cols <- c(sample_names, input_names)
readDepth <- apply(counts[, (cols), with= F], 2, sum)
normPseudo <- readDepth/min(readDepth)*pseudocount
normPseudo <- round(normPseudo, 3)
print(paste("Normalized pseudocounts:", paste0(normPseudo, collapse = ", ")))
counts[, (cols):= lapply(seq(.SD), function(i) ifelse(.SD[[i]]==0, .SD[[i]]+normPseudo[i], .SD[[i]])) , .SDcols= cols]

# Save filtered counts file ----
fwrite(counts,
       filtered_output,
       col.names = TRUE,
       row.names = FALSE,
       sep = "\t",
       quote= F,
       na = NA)
