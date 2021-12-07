#' Process ChIP-Seq/ATAC-Seq data
#'
#' Starts from fastq files, outputs bam, bam statistics, bed containing collapsed reads and bw
#'
#' @param fq1 fastq file path (read1)
#' @param fq2 fastq file path (read2). Default is null -> fq1 treated as single end!
#' @param chrom_sizes data.table containing "seqnames" and "seqlengths" columns. 
#' @param basename_prefix prefix to add at the beginning of all output files
#' @param Rsubread_index_prefix Rsubread index to use for alignment. Attached seqnames and seqlenths will be used.
#' @param bam_folder Folder for output bam files.
#' @param bw_folder Folder for output bw files.
#' @param extend length to which single end reads should be extended (bed and bw)
#' @param maxMismatches default 2
#' @param nTrim3 default 0
#' @param nTrim5 default 0
#' @param use_samtools If FALSE, avoids using samtools by saving alignments in a non-binary sam file. default= TRUE.
#' @examples require(data.table)
#' require(GenomicRanges)
#' require(Rsubread)
#' vl_ChIP_pipeline(fq1= "/groups/stark/vloubiere/projects/available_data_dm3/test/SRR1042411.fastq.gz",
#' basename_prefix = "test_single_end",
#' Rsubread_index_prefix = "/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/dm3_genome/index/dm3_genome_idx",
#' bam_folder = "/groups/stark/vloubiere/projects/available_data_dm3/test/",
#' bw_folder = "/groups/stark/vloubiere/projects/available_data_dm3/test/")
#' vl_ChIP_pipeline(fq1= "/groups/stark/vloubiere/projects/available_data_dm3/db/fastq/SRR7813063_1.fastq.gz",
#' fq2= "/groups/stark/vloubiere/projects/available_data_dm3/db/fastq/SRR7813063_2.fastq.gz",
#' basename_prefix = "test_paired_end",
#' Rsubread_index_prefix = "/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/dm3_genome/index/dm3_genome_idx",
#' bam_folder = "/groups/stark/vloubiere/projects/available_data_dm3/test/",
#' bw_folder = "/groups/stark/vloubiere/projects/available_data_dm3/test/")
#' @return bam , bam_stats.txt, uniq.bed and .bw files
#' @export


vl_ChIP_pipeline <- function(fq1,
                             fq2= NULL,
                             chrom_sizes,
                             basename_prefix= "",
                             Rsubread_index_prefix= NULL,
                             bam_folder= NULL,
                             bw_folder,
                             extend= 300,
                             max_frag_size_pe= 500,
                             maxMismatches= 2,
                             nTrim3= 0,
                             nTrim5= 0)
{
  if(!is.data.table(chrom_sizes))
    stop("chrom sizes must be a data.table object")
  if(!all(c("seqnames", "seqlengths") %in% colnames(chrom_sizes)))
     stop("chrom sizes data.table object should contain seqnames and seqlengths columns!")
  if(any(!grepl(".fq$|.fastq$|.fq.gz$|.fastq.gz$", fq1)))
    stop("Unsupported fq extension. Supported: fq|fastq|fq.gz|fastq.gz")
     
  # Compute output name
  .bn <- gsub(".fq$|.fastq$|.fq.gz$|.fastq.gz$", "", basename(fq1))
  if(!is.null(fq2))
    .bn <- gsub("_1$", "", .bn)
  output <- paste0(bam_folder, basename_prefix, .bn, ".sam")
  # Alignment
  print("START alignment!")
  Rsubread::align(index = Rsubread_index_prefix,
                  readfile1 = fq1, 
                  readfile2 = fq2, 
                  type= "dna", 
                  output_file = output, 
                  maxMismatches= maxMismatches, 
                  unique= T, 
                  nTrim3 = nTrim3,
                  nTrim5 = nTrim5, 
                  nthreads= getDTthreads()-2, 
                  output_format = "SAM")
  
  # Import reads
  print("Import reads!")
  reads <- data.table::fread(output, 
                             fill= T, 
                             select = c("V1", "V2", "V3", "V4", "V5", "V10"), 
                             col.names = c("ID", "flag", "seqnames", "read_most_left_pos", "mapq", "read"))
  reads <- reads[seqnames %in% chrom_sizes$seqnames]
  reads[, read_length:= nchar(read)]
  reads$read <- NULL
  cols <- c("flag", "read_most_left_pos", "mapq", "read_length")
  reads[, (cols):= lapply(.SD, as.numeric), .SDcols= cols]

  # Clean
  print("Import done -> start cleaning reads!")
  reads <- reads[mapq>=30]
  if(is.null(fq2))
  {
    reads[flag==0, start:= read_most_left_pos]
    reads[flag==16, start:= read_most_left_pos+read_length-1-extend]
    reads <- unique(reads[, .(seqnames, start)])
    reads[, end:= start+extend]
    reads[start<1, start:= 1]
    reads[chrom_sizes, end:= i.seqlengths, on= c("seqnames", "end>seqlengths")]
    res <- GenomicRanges::GRanges(reads)
  }else
  {
    res <- unique(reads[flag %in% c(99, 163, 83, 147), .(ID, seqnames)])
    res[reads[flag %in% c(99, 163)], start:= read_most_left_pos, on= "ID"]
    res <- na.omit(res)
    res[reads[flag %in% c(83, 147)], end:= read_most_left_pos+read_length-1, on= "ID"]
    res <- unique(na.omit(res[(end-start)<max_frag_size_pe, !"ID"]))
    res <- GenomicRanges::GRanges(res)
  }
  # Export result as bed
  print("Export bed!")
  rtracklayer::export.bed(res, 
                          gsub(".bam$|.sam$", "_uniq.bed", output))
  # Generate and export bw
  print("Generate bw!")
  total_reads <- length(res)
  cov <- GenomicRanges::coverage(res)/total_reads*1e6
  bw_output <- paste0(bw_folder, basename_prefix, .bn, ".bw")
  rtracklayer::export.bw(GRanges(cov), con= bw_output)
  
  # gzip sam file
  system(paste("gzip", output))
}

