#' ORFtag pipeline
#' 
#' Takes as input a (corerectly formated) metadata file, save the processed metadata file and returns the command lines to 1/ extract reads from VBC bam file, 2/ trim the reads and align to mouse/human genome (see 'species' column of the metadata table) and return alignment statistics as well as collapsed reads and 4/ assign insertions to closest downstream genes.
#'
#' @param metadata The path to a .txt, tab-separated metadata file containing at least 12 columns. See vlfunctions::vl_metadata_ORFtag for an example.
#' @param processed_metadata_output An .rds path where to save the processed metadata file (containing the paths of output files).
#' @param Rpath Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript" 
#' @param cores Number of cores per job. Default= 8
#' @param mem Memory per job (in Go). Default= 32.
#' @param overwrite Should existing files be overwritten?
#' @param submit Should the command be submitted? default= G
#' @param wdir The working directory to use. defaut= getwd().
#' @param logs Path to save logs. Default= "db/logs"
#'
#' @return Return a data.table containing, for each sampleID, the concatenated commands that are required to process the data. These commands can then be submitted using ?vl_bsub().
#' @export
#'
#' @examples
#' A metadata example can be found at ?vl_metadata_ORFtag

vl_ORFtag_pipeline <- function(metadata,
                               processed_metadata_output,
                               Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                               cores= 8,
                               mem= 32,
                               overwrite= F,
                               submit= F,
                               wdir= getwd(),
                               logs= "db/logs/ORFtag")
{
  # Import metadata and check format ----
  meta <- fread(metadata, header = T)
  cols <- c("user", "batch", "screen", "condition", "replicate", "barcodes", "sampleID", "expName", "species", "layout", "sequencer", "bam_path")
  if(!all(cols %in% names(meta)))
    stop(paste(c("metadata file should contain at least these 12 columns:", cols), collapse = " "))
  if(!all(meta[, sampleID==paste0(screen, "_", condition, "_rep", replicate)]))
    stop("sampleID should be the catenation of screen, condtion and replicate. Sany samplex with the same sampleID will be collapsed!")
  
  # Generate output paths ----
  meta[, fq1:= paste0("db/fastq/ORFtag/", gsub(".bam", "", basename(bam_path)), "_", make.unique(sampleID), "_1.fq.gz")]
  meta[(layout=="PAIRED"), fq2:= gsub("_1.fq.gz$", "_2.fq.gz", fq1)]
  meta[, fq1_trim:= gsub(".fq.gz$", "_trimmed.fq.gz", fq1)]
  meta[, bam:= paste0("db/bam/ORFtag/", sampleID, ".bam")] # re-sequencing are merged from this step on!
  meta[, bam_stats:= gsub(".bam$", "_stats.txt", bam)]
  meta[, bam_unique:= paste0("db/bam_unique/ORFtag/", sampleID, "_q30_unique.bam")]
  meta[, bed_file:= paste0("db/bed/ORFtag/", sampleID, ".bed")]
  meta[, counts_same_strand:= paste0("db/gene_assignment/ORFtag/", sampleID, "_same_strand.txt")]
  meta[, counts_rev_strand:= paste0("db/gene_assignment/ORFtag/", sampleID, "_rev_strand.txt")]
  
  # Save processed metadata ----
  if(!grepl(".rds$", processed_metadata_output))
    stop("processed_metadata_output should end up with a .rds extension") else
      saveRDS(meta, processed_metadata_output)
  
  # Create output directories ----
  dirs <- c(logs,
            na.omit(unique(dirname(unlist(meta[, fq1:counts_rev_strand])))))
  if(any(!dir.exists(dirs)))
  {
    dirs <- dirs[!dir.exists(dirs)]
    # Ask the user for their name
    ask <- readline(prompt = paste0(c("Create output directories: ", dirs, "? (yes/cancel) "), collapse = " "))
    if(ask=="yes")
      sapply(dirs, dir.create, showWarnings = F, recursive = T) else
        stop("Stopped because output directories were missing.")
  }
  
  # Load modules ----
  meta[, load_cmd:= paste(c(paste("cd", wdir),
                            "module load build-env/2020",
                            "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                            "module load samtools/1.9-foss-2018b",
                            "module load bowtie2/2.3.4.2-foss-2018b"), collapse = "; ")]
  
  # Demultiplex VBC bam file ----
  meta[, demultiplex_cmd:= {
    if(.N!=1)
      stop("Unique bam file should be provided!")
    if(overwrite | !file.exists(fq1)) # fq2 is not used anyway
    {
      BC <- paste0("'^BC:Z:", paste0(unlist(tstrsplit(barcodes, "\\|")), collapse = "|^BC:Z:"), "'")
      fq_prefix <- gsub("_1.fq.gz$", "", fq1)
      cmd <- paste("samtools view -@", cores-1, bam_path, 
                   "| perl", 
                   fcase(layout=="SINGLE" & sequencer=="NextSeq", # se reads, BC is in column 12
                         system.file("ORFtag_pipeline", "demultiplex_se_12.pl", package = "vlfunctions"), 
                         layout=="PAIRED" & sequencer=="NextSeq", # pe reads, BC is in column 12 (typically what we use)
                         system.file("ORFtag_pipeline", "demultiplex_pe_12.pl", package = "vlfunctions"), 
                         layout=="SINGLE" & sequencer=="HiSeq", # pe reads, BC is in column 12
                         system.file("ORFtag_pipeline", "demultiplex_se_14.pl", package = "vlfunctions")),
                   BC,
                   fq_prefix,
                   "; gzip",
                   paste0(fq_prefix, "_1.fq"))
      if(!is.na(fq2))
        cmd <- paste0(cmd, "; gzip ", fq_prefix, "_2.fq")
      cmd
    }
  }, .(layout, barcodes, fq1, fq2)]
  
  # Trim reads ----
  meta[, trim_cmd:= {
    if(overwrite | !file.exists(fq1_trim))
      paste0("trim_galore --gzip -o ", dirname(fq1_trim), "/ ", fq1)
  }, fq1_trim]
  
  # Alignment ----
  meta[, aln_cmd:= {
    if(overwrite | !file.exists(bam))
    {
      #### BOWTIE 2
      x <- switch(species,
                  "mouse" = "/groups/stark/vloubiere/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome",
                  "human" = "/groups/stark/vloubiere/genomes/Homo_sapiens/hg38/Bowtie2Index/genome")
      paste("bowtie2 -p", cores,
            "-U", paste0(fq1_trim, collapse= ","),
            "-x", x,
            "| samtools sort -@", cores-1, "-o", bam)
    }
  }, .(bam, species)]
  
  # Get bam stats ----
  meta[, stats_cmd:= {
    if(overwrite | !file.exists(bam_stats))
      paste("samtools stats -@", cores-1, bam, "| grep ^SN | cut -f 2- >", bam_stats)
  }, .(bam, bam_stats)]
  
  # Collapsed bam files ----
  meta[, collapse_cmd:= {
    if(overwrite | !file.exists(bam_unique))
      paste("samtools sort -n -@", cores-1, bam, 
            "| samtools fixmate -m - - | samtools sort -@", cores-1, 
            "| samtools markdup -r - - | samtools view -q 30 -b -o",  bam_unique)
  }, .(bam, bam_unique)]
  
  # Compute insertions ----
  meta[, insertions_cmd:= {
    if(overwrite | any(!file.exists(c(counts_same_strand, counts_rev_strand, bed_file))))
    {
      x <- switch(species,
                  "mouse" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
                  "human" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_hg38.gtf")
      paste(Rpath,
            system.file("ORFtag_pipeline", "bamToBed_and_assign_insertions.R", package = "vlfunctions"), 
            bam_unique,
            x,
            bed_file,
            gsub("_same_strand.txt$", "", counts_same_strand))
    }
  }, .(species, bam_unique, bed_file, counts_same_strand, counts_rev_strand)]
  
  # Return commands ----
  cols <- c("load_cmd", "demultiplex_cmd", "trim_cmd", "aln_cmd", "stats_cmd", "collapse_cmd", "insertions_cmd")
  cols <- cols[cols %in% names(meta)]
  cmd <- meta[, {
    cmd <- paste0(unique(na.omit(unlist(.SD))), collapse = "; ")
    if(cmd==load_cmd[1])
      cmd <- as.character()
    .(cmd= cmd)
  }, sampleID, .SDcols= cols]
  
  # Submit commands ----
  if(nrow(cmd))
  {
    if(submit)
      cmd[, {
        vl_bsub(cmd, 
                cores= cores, 
                m = mem, 
                name = "ORFtag", 
                t = '1-00:00:00',
                o= normalizePath(logs),
                e= normalizePath(logs))
      }, cmd] else
        return(cmd)
  }else
    print("All output files already existed! No command submitted ;)")
}

#' ORFtrap_call_hits
#' 
#' @description Function to call hits from a given screen. Receives sorted and unsorted counts as input, computes FC table 
#' 
#' @param sorted_forward_counts Sorted forward counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param unsorted_forward_counts Unsorted (input) forward counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param exons_start_gtf gtf exon file. Used for consistent ordering of output FC table
#' @param name Name to be appended at the beginning of output file
#' @param output_suffix Suffix to be appended at the end of output file. Default to "_vs_unsorted.txt"
#' @param output_folder_FC_file Output folder for FC files
#'
#' @return Returns FC tables containing DESeq2-like columns
#' @export
#'
#' @examples
#' vl_ORFtrap_call_hits(sorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_sort_rep1_same_strand.txt",
#'                                                "db/gene_assignment/ORFtag/Activator2_sort_rep2_same_strand.txt"),
#'                      unsorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_input_rep1_same_strand.txt",
#'                                                  "db/gene_assignment/ORFtag/Activator2_input_rep2_same_strand.txt"),
#'                      exons_start_gtf = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
#'                      name = "Activator_2",
#'                      output_suffix = "_vs_input.txt",
#'                      output_folder_FC_file = "db/FC_tables/ORFtag/Activator_2")

vl_ORFtrap_call_hits <- function(sorted_forward_counts, 
                              unsorted_forward_counts,
                              exons_start_gtf,
                              name,
                              output_suffix= "_vs_unsorted.txt",
                              output_folder_FC_file)
{
  require(rtracklayer)
  require(data.table)
  require(GenomicRanges)
  
  # Checks
  if(!is.character(name) | length(name)!=1)
    stop("name should be a unique character specifying the name of the screen")
  if(anyDuplicated(sorted_forward_counts))
    stop("Duplicated filenames in sorted_forward_counts")
  if(anyDuplicated(unsorted_forward_counts))
    stop("Duplicated filenames in unsorted_forward_counts")
  if(!is.character(output_suffix) | length(output_suffix)!=1)
    stop("output_suffix should be a unique character specifying the name of the screen")
  
  # Import counts
  input <- rbindlist(lapply(unsorted_forward_counts, fread))
  sample <- rbindlist(lapply(sorted_forward_counts, fread))
  dat <- rbindlist(list(count.input= input, count.sample= sample), idcol= "cdition")
  dat <- na.omit(dat[dist<2e5])
  total <- dat[, .N, .(cdition)]
  dat <- dat[, .(count= .N), .(gene_id, cdition)]
  dat <- dcast(dat, gene_id~cdition, value.var = "count")
  
  # Import gene exons
  genes <- rtracklayer::import(exons_start_gtf)
  genes <- as.data.table(mcols(genes)[, c("gene_id", "gene_name")])
  genes <- unique(genes)
  setorderv(genes, "gene_id")
  
  # Format dat
  dat <- merge(genes, dat, by= "gene_id", all.x= T)
  dat[is.na(count.input), count.input:= 0]
  dat[is.na(count.sample), count.sample:= 0]
  dat[, total.input:= total[cdition=="count.input", N]]
  dat[, total.sample:= total[cdition=="count.sample", N]]
  
  # Fisher test sample vs input
  dat[count.sample>=3, c("OR", "pval"):= {
    .t <- matrix(c(total.input-count.input+1,
                   count.input+1,
                   total.sample-count.sample+1,
                   count.sample+1),
                 ncol= 2)
    fisher.test(.t, alternative = "greater")[c("estimate", "p.value")]
  }, .(count.input, total.input, count.sample, total.sample)]
  dat[count.sample>=3, log2OR:= log2(OR)]
  dat[count.sample>=3, padj:= p.adjust(pval, "fdr")]
  dat[, hit:= padj<0.001 & log2OR>=1]
  
  # Clean and save
  dat <- genes[dat, on= c("gene_id", "gene_name")]
  dat$OR <- dat$pval <- NULL
  dir.create(output_folder_FC_file, showWarnings = F, recursive = T)
  FC_table <- paste0(output_folder_FC_file, "/", name, output_suffix)
  fwrite(dat, 
         FC_table, 
         col.names = T, 
         row.names = F, 
         sep= "\t",
         quote= F, 
         na= NA)
  return(paste0(name, ": ", sum(dat$hit, na.rm = T), " hits were called!\nFC file -> ", FC_table, "\n"))
}

#' vl_ORFtrap_call_hits_strandBias
#' 
#' @description Function to call hits from a given screen. Receives sorted and unsorted counts as input, computes FC table 
#' 
#' @param sorted_forward_counts Sorted forward counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param sorted_reverse_counts Sorted reverse counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param unsorted_forward_counts Unsorted (input) forward counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param unsorted_reverse_counts Unsorted (input) reverse counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param exons_start_gtf gtf exon file. Used for consistent ordering of output FC table
#' @param name Name to be appended at the beginning of output file
#' @param output_suffix Suffix to be appended at the end of output file. Default to "_vs_revStrand.txt"
#' @param output_folder_FC_file Output folder for FC files
#'
#' @return Returns FC tables containing DESeq2-like columns
#' @export
#'
#' @examples
#' vl_ORFtrap_call_hits_strandBias(sorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_sort_rep1_same_strand.txt",
#'                                                           "db/gene_assignment/ORFtag/Activator2_sort_rep2_same_strand.txt"),
#'                                 sorted_reverse_counts = c("db/gene_assignment/ORFtag/Activator2_sort_rep1_rev_strand.txt",
#'                                                           "db/gene_assignment/ORFtag/Activator2_sort_rep2_rev_strand.txt"),
#'                                 unsorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_input_rep1_same_strand.txt",
#'                                                             "db/gene_assignment/ORFtag/Activator2_input_rep2_same_strand.txt"),
#'                                 unsorted_reverse_counts = c("db/gene_assignment/ORFtag/Activator2_input_rep1_rev_strand.txt",
#'                                                             "db/gene_assignment/ORFtag/Activator2_input_rep2_rev_strand.txt"),
#'                                 exons_start_gtf = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
#'                                 name = "Activator_2",
#'                                 output_suffix = "_vs_input.txt",
#'                                 output_folder_FC_file = "db/FC_tables/ORFtag/")

vl_ORFtrap_call_hits_strandBias <- function(sorted_forward_counts, 
                                         sorted_reverse_counts, 
                                         unsorted_forward_counts,
                                         unsorted_reverse_counts,
                                         exons_start_gtf,
                                         name,
                                         output_suffix= "_vs_revStrand",
                                         output_folder_FC_file)
{
  require(rtracklayer)
  require(data.table)
  require(GenomicRanges)
  
  # Checks
  if(!is.character(name) | length(name)!=1)
    stop("name should be a unique character specifying the name of the screen")
  if(anyDuplicated(sorted_forward_counts))
    stop("Duplicated filenames in sorted_forward_counts")
  if(anyDuplicated(unsorted_forward_counts))
    stop("Duplicated filenames in unsorted_forward_counts")
  if(!is.character(output_suffix) | length(output_suffix)!=1)
    stop("output_suffix should be a unique character specifying the name of the screen")
  
  # Import counts
  inputFw <- rbindlist(lapply(unsorted_forward_counts, fread))
  inputRev <- rbindlist(lapply(unsorted_reverse_counts, fread))
  sampleFw <- rbindlist(lapply(sorted_forward_counts, fread))
  sampleRev <- rbindlist(lapply(sorted_reverse_counts, fread))
  
  dat <- rbindlist(list(count.sample.fw= sampleFw,
                        count.sample.rev= sampleRev,
                        count.input.fw= inputFw,
                        count.input.rev= inputRev),
                   idcol= "cdition")
  dat <- na.omit(dat[dist<5e4])
  dat <- dat[, .(count= .N), .(gene_id, cdition)]
  dat <- dcast(dat,
               gene_id~cdition,
               value.var = "count",
               fill= 0)
  
  # Import gene exons
  genes <- rtracklayer::import(exons_start_gtf)
  genes <- as.data.table(mcols(genes)[, c("gene_id", "gene_name")])
  genes <- unique(genes)
  setorderv(genes, "gene_id")
  
  # Format dat
  dat <- merge(genes, dat, by= "gene_id", all.x= T)
  
  # Fisher test fw vs rev
  dat[count.sample.fw>=3, c("OR", "pval"):= {
    .t <- matrix(c(count.sample.fw+1,
                   count.sample.rev+1,
                   count.input.fw+1,
                   count.input.rev+1),
                 ncol= 2)
    fisher.test(.t, alternative = "greater")[c("estimate", "p.value")]
  }, .(count.sample.fw, count.sample.rev, count.input.fw, count.input.rev)]
  dat[count.sample.fw>=3, log2OR:= log2(OR)]
  dat[count.sample.fw>=3, padj:= p.adjust(pval, "fdr")]
  dat[, hit:= padj<0.05 & log2OR>=1]
  
  # Add binomial test
  dat[count.sample.fw>=3, binom_pval:= {
    binom.test(x = count.sample.fw+1,
               n = count.sample.fw+count.sample.rev+1,
               p = (count.input.fw+1)/(count.input.fw+count.input.rev+1),
               alternative = "greater")["p.value"]
  }, .(count.sample.fw, count.sample.rev, count.input.fw, count.input.rev)]
  dat[count.sample.fw>=3, binom_padj:= p.adjust(binom_pval, "fdr")]
  dat[, binom_hit:= binom_padj<0.001]
  
  # Clean and save
  dat <- genes[dat, on= c("gene_id", "gene_name")]
  dat$OR <- dat$pval <- dat$binom_pval <- NULL
  dir.create(output_folder_FC_file, showWarnings = F, recursive = T)
  FC_table <- paste0(output_folder_FC_file, "/", name, output_suffix)
  fwrite(dat, 
         FC_table, 
         col.names = T, 
         row.names = F, 
         sep= "\t",
         quote= F, 
         na= NA)
  return(paste0(name, ": ", sum(dat$hit, na.rm = T), " hits were called!\nFC file -> ", FC_table, "\n"))
}
