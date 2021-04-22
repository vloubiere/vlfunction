#' Process ChIP-Seq/ATAC-Seq data
#'
#' Starts from fastq files, outputs bam, bam statistics, bed containing collapsed reads and bw
#'
#' @param fq1 fastq file path (read1)
#' @param fq2 fastq file path (read2). Default is null -> fq1 treated as single end!
#' @param basename_prefix prefix to add at the beginning of all output files
#' @param Rsubread_index_prefix Rsubread index to use for alignment. Attached seqnames and seqlenths will be used.
#' @param bam_folder Folder for output bam files.
#' @param bw_folder Folder for output bw files.
#' @param extend length to which single end reads should be extended (bed and bw)
#' @param maxMismatches default 2
#' @param nTrim3 default 0
#' @param nTrim5 default 0
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
                             basename_prefix= "",
                             Rsubread_index_prefix= NULL,
                             bam_folder= NULL,
                             bw_folder,
                             extend= 300,
                             maxMismatches= 2,
                             nTrim3= 0,
                             nTrim5= 0)
{
  # Compute output name
  .bn <- strsplit(basename(fq1), "[.]")[[1]][1]
  output <- paste0(bam_folder, basename_prefix, "_", .bn, ".bam")
  # Alignment
  print("START alignment!")
  stats <- capture.output(Rsubread::align(index = Rsubread_index_prefix,
                                          readfile1 = fq1, 
                                          readfile2 = fq2, 
                                          type= "dna", 
                                          output_file = output, 
                                          maxMismatches= maxMismatches, 
                                          unique= T, 
                                          nTrim3 = nTrim3,
                                          nTrim5 = nTrim5, 
                                          nthreads= getDTthreads()-2))
  writeLines(stats, gsub(".bam$", "_stats.txt", output))
  # Chromosome infos
  chr <- data.table::fread(list.files(dirname(Rsubread_index_prefix), ".reads$", full.names = T), 
                           col.names = c("length", "seqnames"))
  # Import reads
  print("Import reads!")
  cmd <- paste("module load build-env/2020; module load samtools/1.9-foss-2018b; samtools view -@", 
               getDTthreads()-1, 
               output)
  reads <- data.table::fread(cmd= cmd, fill= T)
  # Clean
  print("Import done -> start cleaning reads!")
  reads <- reads[V3 %in% chr$seqnames & V5>30]
  if(is.null(fq2))
  {
    reads <- unique(reads[V2 %in% c(0, 16), 
                          .(V2, V3, V4, V10= nchar(V10))]) # nchar(V10)= read lengths
    left <- reads[V2==0, .(seqnames= V3, 
                           start= V4,
                           end= V4+extend)]
    right <- reads[V2==16, .(seqnames= V3,
                             end= V4+V10)]
    right[, start:= end-300]               
    reads <- rbind(left, right)
    reads[start<1, start:= 1]
    reads[chr, end:= i.length, on= c("seqnames", "end>length")]
  }else
  {
    reads <- reads[V2 %in% c(83, 99, 147, 163), 
                   .(V1, V2, V3, V4, V10= nchar(V10))] # nchar(V10)= read lengths
    data.table::setkeyv(reads, "V1")
    reads <- reads[, .(seqnames= V3, 
                       start= min(V4), 
                       end= max(V4)+V10), V1]
    reads <- unique(reads[, seqnames:end])
  }
  # Export result as bed
  print("Export bed!")
  reads <- GRanges(reads)
  rtracklayer::export.bed(reads, 
                          gsub(".bam$", "_uniq.bed", output))
  # Generate and export bw
  print("Generate bw!")
  total_reads <- length(reads)
  cov <- GenomicRanges::coverage(reads)/total_reads*1e6
  output <- paste0(bw_folder, basename_prefix, "_", .bn, ".bw")
  rtracklayer::export.bw(GRanges(cov), con= output)
}

