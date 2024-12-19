#' ORFtag pipeline
#' 
#' To use this this pipeline, please install the vl_function package using install_github("vloubiere/vlfunction")
#' or devtools::install("/groups/stark/vloubiere/vlfunction/") in your interactive Rstudio session AND
#' in the local installation of R 4.3.0 ('/software/f2022/software/r/4.3.0-foss-2022b/bin/R') 
#' which will be used to submit R sub-script.
#' 
#' Takes as input a (correctly formatted) metadata file, saves the processed metadata file and returns and.or submit the command lines to: \cr
#'  1/ extract reads from VBC bam file \cr
#'  2/ trim the reads \cr
#'  3/ Aligns to mouse/human genome (see 'genome' column of the metadata table) and returns a bam file \cr
#'  4/ Return alignment statistics \cr
#'  5/ Collapse unique reads and stores them into a bam file \cr
#'  6/ assign insertions to closest downstream genes \cr
#'  
#'  Then, the second function vl_ORFtrap_call_hits() can be used to call confident hits.
#'
#' @param metadata The path to a correctly formatted .xlsx, .rds or .txt metadata file, or a data.table. See examples for a template.
#' @param processed_metadata_output A .rds path to save the processed metadata file, containing the output file paths. By default, "_processed.rds" will be appended to the excel file path. 
#' @param fq_output_folder Output folder for demultiplexed fq files. Default= "/scratch/stark/vloubiere/fq/ORFtag/".
#' @param bam_output_folder Output folder for aligned bam files Default= "/scratch/stark/vloubiere/bam/ORFtag/".
#' @param alignment_stats_folder Output folder for alignment statistics. Default= "db/alignment_stats/ORFtag/".
#' @param bam_unique_output_folder Output folder for .bam files containing unique insertions. Default= "db/bam_unique/ORFtag/".
#' @param bed_output_folder Output folder for .bed files containing unique insertions. Default= "db/bed/ORFtag/".
#' @param exon_assignment_output_folder Output folder for counts_same_strand and counts_rev_strand files, containing exon assignment for each unique insertion. Default= "db/exon_assignment/ORFtag/".
#' @param tmp_folder Output folder for temporary bam files during sorting. Default= "/scratch/stark/vloubiere/ORFtag/tmp/sort/".
#' @param Rpath Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of cores per job. Default= 8.
#' @param mem Memory per job (in Go). Default= 32.
#' @param overwrite Should existing files be overwritten? Default= FALSE, meaning that only the commands for which output files do not exist will be returned and/or submitted.
#' @param submit Should the command be submitted? Default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs".
#' @param time The time required for the SLURM scheduler. Default= '1-00:00:00'.
#'
#' @return Return a data.table containing, for each sampleID, the concatenated commands that are required to process the data.
#' These commands can then be submitted either directly via the function, or using vl_bsub()....
#'
#' @examples
#' library(vlfunctions)
#' 
#' # Example metadata ----
#' metadata <- system.file("ORFtag_pipeline", "metadata_ORFtag.xlsx", package = "vlfunctions")
#' 
#' # Generate command lines ----
#' cmd <- vl_ORFtag_pipeline(metadata= metadata,
#'                           overwrite= TRUE,
#'                           submit= FALSE)
#'                           
#' # To actually submit the jobs, switch to submit= TRUE ----
#' vl_ORFtag_pipeline(metadata= metadata,
#'                    processed_metadata_output = "Rdata/metadata_ORFtag_processed.rds",
#'                    overwrite= FALSE, # Ensures that existing files are not processed again 
#'                    submit= TRUE)
#'                    
#' # Processed metadata
#' processed <- readRDS("Rdata/metadata_ORFtag_processed.rds")
#' 
#' # To call hits
#' see ?vl_ORFtrap_call_hits()
#' 
#' # To create .gtf files containing non-first exon start coordinates ------------------
#' # Mouse (mm10):
#' # gtf <- rtracklayer::import("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz")
#' gtf <- rtracklayer::import("db/gtf/gencode.vM25.basic.annotation.gtf.gz")
#' exons <- gtf[gtf$transcript_type=="protein_coding" & gtf$type=="exon" & gtf$exon_number>1]
#' exons <- GenomicRanges::resize(exons, 1, "start")
#' mcols(exons) <- mcols(exons[, c("gene_id", "gene_name", "mgi_id", "exon_number", "exon_id")])
#' exons$gene_id <- gsub("[.*]..*", "\\1", exons$gene_id)
#' exons$exon_id <- gsub("[.*]..*", "\\1", exons$exon_id)
#' exons <- unique(exons)
#' rtracklayer::export(exons,
#'                     "db/gtf/exons_start_mm10.gtf")
#'                     
#' # Human (hg38):
#' # gtf <- rtracklayer::import("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz")
#' 
#' @export
vl_ORFtag_pipeline <- function(metadata,
                               processed_metadata_output,
                               fq_output_folder= "/scratch/stark/vloubiere/fq/ORFtag/",
                               bam_output_folder= "/scratch/stark/vloubiere/bam/ORFtag/",
                               alignment_stats_folder= "db/alignment_stats/ORFtag/",
                               bam_unique_output_folder= "db/bam_unique/ORFtag/",
                               bed_output_folder= "db/bed/ORFtag/",
                               exon_assignment_output_folder= "db/exon_assignment/ORFtag/",
                               tmp_folder= "/scratch/stark/vloubiere/ORFtag/tmp/sort/",
                               Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                               cores= 8,
                               mem= 32,
                               overwrite= FALSE,
                               submit= FALSE,
                               wdir= getwd(),
                               logs= "db/logs/ORFtag",
                               time= '1-00:00:00')
{
  # Import metadata and check layout ----
  meta <- if(is.data.table(metadata))
    data.table::copy(metadata) else
      vl_import_metadata_sheet(metadata)
  
  # Check required columns ----
  # Design
  missing_cols <- setdiff(c("screen", "sampleID", "genome", "layout", "sequencer", "fq1", "fq2"), names(meta))
  if(length(missing_cols))
    stop(paste("Metadata required columns mising:", paste0(missing_cols, collapse = ", ")))
  # Check fq files
  meta[, fq1:= as.character(fq1)]
  meta[, fq2:= as.character(fq2)]
  if(nrow(meta[is.na(fq1)]) | nrow(meta[layout=="PAIRED" & is.na(fq2)]))
  {
    missing_cols <- setdiff(c("bam_path", "barcodes"), names(meta))
    if(length(missing_cols))
      stop(paste("fq1 files not provided and the metadata is missing the following columns for demultiplexing:",
                 paste0(missing_cols, collapse = ", ")))
  }
  
  # Check arguments ----
  # If fq files provided, make sure they exist
  fqs <- na.omit(unlist(meta[!is.na(fq1)|!is.na(fq2), .(fq1, fq2)]))
  if(length(fqs))
  {
    if(any(!grepl(".fq.gz$", fqs)))
      stop("Some user-provided .fq files do not have the correct extesion, .fq.gz")
    if(any(!file.exists(fqs)))
      stop("Some user-provided .fq files do not exist. Check that the path is correct")
  }
  # Check metadata output
  if(submit && !grepl(".rds$", processed_metadata_output))
    stop("processed_metadata_output should end up with a .rds extension")
  # Alignment indexes
  if(any(!meta$genome %in% c("mm10", "hg38")))
    stop("Only mm10 and hg38 are supported! For other genomes, please provide path to the corresponding bowtie 2 index.")
  # layout
  if(!all(meta$layout %in% c("SINGLE", "PAIRED")))
    stop("layout column should only contain 'SINGLE' or 'PAIRED'")
  # Warning for duplicated sample IDs
  potential_dup <- meta[, .N, sampleID][N>1]
  if(nrow(potential_dup))
    warning(paste("The following sampleIDs were present more than once and will be merged (resequenced)?:",
                  paste0(potential_dup$sampleID, collapse = ", ")))
  
  # Generate output paths ----
  # Temp files
  meta[is.na(fq1) & !is.na(bam_path), fq1:= {
    paste0(fq_output_folder, "/",
           gsub(".bam", paste0("_", gsub("|", "_", barcodes, fixed = T), "_1.fq.gz"), basename(bam_path)))
  }, .(barcodes, bam_path)]
  meta[is.na(fq2) & layout=="PAIRED", fq2:= gsub("_1.fq.gz$", "_2.fq.gz", fq1)]
  meta[, fq1_trim:= paste0(fq_output_folder, "/", gsub(".fq.gz$", "_trimmed.fq.gz", basename(fq1)))]
  # re-sequenced samples merged from this step on!
  meta[, bam:= paste0(bam_output_folder, "/", sampleID, ".bam")] 
  # Temp files
  meta[, bam_stats:= paste0(alignment_stats_folder, gsub(".bam$", "_stats.txt", basename(bam)))]
  meta[, bam_unique:= paste0(bam_unique_output_folder, sampleID, "_q30_unique.bam")]
  meta[, bed_file:= paste0(bed_output_folder, sampleID, ".bed")]
  meta[, counts_same_strand:= paste0(exon_assignment_output_folder, sampleID, "_same_strand.txt")]
  meta[, counts_rev_strand:= paste0(exon_assignment_output_folder, sampleID, "_rev_strand.txt")]
  
  # Print conditions ----
  cat(paste(length(unique(meta$screen)), "screen(s) detected:\n"))
  meta[, cat(paste0(screen, ":\n\t", paste0(unique(sampleID), collapse = "\n\t"), "\n")), screen]
  
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
                   fcase(layout=="SINGLE" & sequencer=="HiSeq", # pe reads, BC is in column 14 (HiSeq only)
                         system.file("ORFtag_pipeline", "demultiplex_se_14.pl", package = "vlfunctions"),
                         layout=="SINGLE", # se reads, BC is in column 12
                         system.file("ORFtag_pipeline", "demultiplex_se_12.pl", package = "vlfunctions"), 
                         layout=="PAIRED", # pe reads, BC is in column 12 (typically what we use)
                         system.file("ORFtag_pipeline", "demultiplex_pe_12.pl", package = "vlfunctions")),
                   BC,
                   fq_prefix,
                   "; gzip -f", paste0(fq_prefix, "_1.fq"))
      if(!is.na(fq2))
        cmd <- paste0(cmd, "; gzip -f ", fq_prefix, "_2.fq")
      cmd
    }
  }, .(layout, fq1, fq2)]
  
  # Trim reads ----
  meta[, trim_cmd:= {
    if(overwrite | !file.exists(fq1_trim))
      paste0("trim_galore --gzip -o ", dirname(fq1_trim), "/ ", fq1)
  }, fq1_trim]
  
  # Alignment ----
  meta[, aln_cmd:= {
    if(overwrite | !file.exists(bam))
    {
      # BOWTIE 2
      x <- switch(genome,
                  "mm10" = "/groups/stark/vloubiere/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome",
                  "hg38" = "/groups/stark/vloubiere/genomes/Homo_sapiens/hg38/Bowtie2Index/genome")
      paste("bowtie2 -p", cores,
            "-U", paste0(fq1_trim, collapse= ","),
            "-x", x,
            "| samtools sort -@", cores-1, "-o", bam)
    }
  }, .(bam, genome)]
  
  # Get bam stats ----
  meta[, stats_cmd:= {
    if(overwrite | !file.exists(bam_stats))
      paste("samtools stats -@", cores-1, bam, "| grep ^SN | cut -f 2- >", bam_stats)
  }, .(bam, bam_stats)]
  
  # Collapsed bam files ----
  meta[, collapse_cmd:= {
    if(overwrite | !file.exists(bam_unique))
    {
      paste("samtools sort -n -@", cores-1, "-T", tmp_folder, bam,
            "| samtools fixmate -m - - | samtools sort -@", cores-1, "-T", tmp_folder,
            "| samtools markdup -r - - | samtools view -q 30 -b -o",  bam_unique)
    }
  }, .(bam, bam_unique)]
  
  # Compute insertions ----
  meta[, insertions_cmd:= {
    if(overwrite | any(!file.exists(c(counts_same_strand, counts_rev_strand, bed_file))))
    {
      # Assignment
      gtf <- switch(genome,
                    "mm10" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
                    "hg38" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_hg38.gtf")
      paste(Rpath,
            system.file("ORFtag_pipeline", "bamToBed_and_assign_insertions.R", package = "vlfunctions"), 
            bam_unique,
            gtf,
            bed_file,
            gsub("_same_strand.txt$", "", counts_same_strand))
    }
  }, .(bam_unique, genome, bed_file, counts_same_strand, counts_rev_strand)]
  
  # Return commands ----
  load_cmd <- paste(c("module load build-env/2020",
                      "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                      "module load samtools/1.9-foss-2018b",
                      "module load bowtie2/2.3.4.2-foss-2018b"),
                    collapse = "; ")
  cols <- intersect(c("demultiplex_cmd", "trim_cmd", "aln_cmd", "stats_cmd", "collapse_cmd", "insertions_cmd"),
                    names(meta))
  cmd <- if(length(cols))
    meta[, .(cmd= paste0(c(load_cmd, unique(na.omit(unlist(.SD)))),
                         collapse = "; ")), sampleID, .SDcols= cols][cmd!=load_cmd] else
                           data.table()
  
  # If commands were generated ----
  if(nrow(cmd))
  {
    # If commands will be submitted ----
    if(submit)
    {
      # Save processed metadata ----
      saveRDS(meta[, -c(cols), with= FALSE],
              processed_metadata_output)
      
      # Create output directories ----
      dirs <- c(logs,
                tmp_folder,
                na.omit(unique(dirname(unlist(meta[, fq1:counts_rev_strand])))))
      if(any(!dir.exists(dirs)))
      {
        sapply(dirs, dir.create, showWarnings = F, recursive = T)
        outDirs <- paste0(c(fq_output_folder,
                            bam_output_folder,
                            alignment_stats_folder,
                            bam_unique_output_folder,
                            bed_output_folder,
                            exon_assignment_output_folder), collapse= ", ")
        print(paste("Output directories were created in:", outDirs, "!"))
      }
      
      # Submit commands ----
      cmd[, {
        vl_bsub(cmd, 
                cores= cores, 
                m = mem, 
                name = "ORFtag", 
                t = time,
                o= paste0(normalizePath(logs), "/", sampleID),
                e= paste0(normalizePath(logs), "/", sampleID),
                wdir = wdir)
      }, .(sampleID, cmd)]
    }else
    {
      return(cmd)
    }
  }else
    warning("All output files already existed! No command submitted ;). Consider overwrite= T if convenient.")
}

#' ORFtrap_call_hits
#' 
#' @description Function to call hits from a given screen. Receives sorted and unsorted counts as input, computes FC table 
#' 
#' @param sorted.forward.counts Character vector of sorted (FACSed) count file paths generated by vl_ORFtag_pipeline().
#' Counts from all provided files will be merged.
#' @param unsorted.forward.counts Character vector of unsorted (input) count file paths generated by vl_ORFtag_pipeline().
#' Counts from all provided files will be merged.
#' @param genome Genome annotation to be used.
#' For now, only "mm10" and "hg38" are supported, and the corresponding non-first-exon .gtf will be used for sorting.
#' See at the bottom of the ?vl_ORFtag_pipeline help page to see how to generate custom non-first-exon start .gtf files for other genomes.
#' @param name Name prefix appended to output file names.
#' @param padj.cutoff The p.adjust cutoff to be used to call hits (>=). Default= 0.001
#' @param log2OR.cutoff The log2OR cutoff to be used to call hits (>=). Default= 1
#' @param log2OR.pseudocount The pseudocount to be added to avoid infinite values. Note that only the log2OR will be affected, not the fisher p.values. Default= 1
#' @param output.suffix Name suffix appended to output file names. Default= "_vs_unsorted.txt".
#' @param output.folder.FC.file Output folder for FC files. Default= "db/FC_tables/ORFtag".
#'
#' @return Returns FC tables containing DESeq2-like columns.
#'
#' @examples
#' library(vlfunctions)
#' 
#' # Call hits ----
#' # 124 hits should be identified with this dataset
#' vl_ORFtrap_call_hits(sorted.forward.counts = c("db/exon_assignment/ORFtag/Activator2_sort_rep1_same_strand.txt",
#'                                                "db/exon_assignment/ORFtag/Activator2_sort_rep2_same_strand.txt"),
#'                      unsorted.forward.counts = c("db/exon_assignment/ORFtag/Activator2_input_rep1_same_strand.txt",
#'                                                  "db/exon_assignment/ORFtag/Activator2_input_rep2_same_strand.txt"),
#'                      genome = "mm10",
#'                      name = "Activator_2",
#'                      output.suffix = "_vs_input.txt")
#' 
#' # Call hits using revese strand (sanity check -> be cautious with the hits that are also found here!) ----
#' # 81 hits should be identified with the reverse strand
#' vl_ORFtrap_call_hits(sorted.forward.counts = c("db/exon_assignment/ORFtag/Activator2_sort_rep1_rev_strand.txt",
#'                                                "db/exon_assignment/ORFtag/Activator2_sort_rep2_rev_strand.txt"),
#'                      unsorted.forward.counts = c("db/exon_assignment/ORFtag/Activator2_input_rep1_rev_strand.txt",
#'                                                  "db/exon_assignment/ORFtag/Activator2_input_rep2_rev_strand.txt"),
#'                      genome = "mm10",
#'                      name = "Activator_2",
#'                      output.suffix = "_vs_input_rev_strand.txt")
#' 
#' # 4 of the hits should be identified with the reverse strand and should be considered carefully
#' hits <- fread("db/FC_tables/ORFtag/Activator_2_vs_input.txt")[(hit), gene_name]
#' sanityCheck <- fread("db/FC_tables/ORFtag/Activator_2_vs_input_rev_strand.txt")[(hit), gene_name]
#' 
#' # Hits to be considered carefully:
#' intersect(hits, sanityCheck)
#' 
#' # To call hits using strand bias (not used)
#' see ?vl_ORFtrap_call_hits_strandBias()
#'                      
#' @export
vl_ORFtrap_call_hits <- function(sorted.forward.counts, 
                                 unsorted.forward.counts,
                                 genome,
                                 name,
                                 padj.cutoff= 0.001,
                                 log2OR.cutoff= 1,
                                 log2OR.pseudocount= 1,
                                 output.suffix= "_vs_unsorted.txt",
                                 output.folder.FC.file= "db/FC_tables/ORFtag/")
{
  require(rtracklayer)
  require(data.table)
  require(GenomicRanges)
  
  # Checks
  if(!is.character(name) | length(name)!=1)
    stop("name should be a unique character specifying the name of the screen")
  if(anyDuplicated(sorted.forward.counts))
    stop("Duplicated filenames in sorted.forward.counts")
  if(anyDuplicated(unsorted.forward.counts))
    stop("Duplicated filenames in unsorted.forward.counts")
  if(!is.character(output.suffix) | length(output.suffix)!=1)
    stop("output.suffix should be a unique character specifying the name of the screen")
  
  # Import counts
  input <- rbindlist(lapply(unsorted.forward.counts, fread))
  sample <- rbindlist(lapply(sorted.forward.counts, fread))
  dat <- rbindlist(list(count.input= input, count.sample= sample), idcol= "cdition")
  dat <- na.omit(dat[dist<2e5])
  total <- dat[, .N, .(cdition)]
  dat <- dat[, .(count= .N), .(gene_id, cdition)]
  dat <- dcast(dat, gene_id~cdition, value.var = "count")
  
  # Import gene exons
  gtf <- switch(genome,
                "mm10" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
                "hg38" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_hg38.gtf")
  genes <- rtracklayer::import(gtf)
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
    .t <- matrix(c(total.input-count.input,
                   count.input,
                   total.sample-count.sample,
                   count.sample),
                 ncol= 2)
    .(fisher.test(.t+log2OR.pseudocount, alternative = "greater")$estimate,
      fisher.test(.t, alternative = "greater")$p.value)
  }, .(count.input, total.input, count.sample, total.sample)]
  dat[count.sample>=3, log2OR:= log2(OR)]
  dat[count.sample>=3, padj:= p.adjust(pval, "fdr")]
  dat[, hit:= padj<=padj.cutoff & log2OR>=log2OR.cutoff]
  
  # Clean and save
  dat <- genes[dat, on= c("gene_id", "gene_name")]
  dat$OR <- dat$pval <- NULL
  dir.create(output.folder.FC.file, showWarnings = F, recursive = T)
  FC_table <- paste0(output.folder.FC.file, "/", name, output.suffix)
  fwrite(dat, 
         FC_table, 
         col.names = T, 
         row.names = F, 
         sep= "\t",
         quote= F, 
         na= NA)
  cat(paste0(name, ": ", sum(dat$hit, na.rm = T), " hits were called!\nFC file -> ", FC_table, "\n"))
}

#' vl_ORFtrap_call_hits_strandBias
#' 
#' @description Function to call hits from a given screen. Receives sorted and unsorted counts as input, computes FC table 
#' 
#' @param sorted.forward.counts Character vector of sorted (FACSed) count file paths generated by vl_ORFtag_pipeline().
#' Counts from all provided files will be merged.
#' @param sorted.reverse.counts Character vector of sorted (FACSed), REVERSED count file paths generated by vl_ORFtag_pipeline().
#' Counts from all provided files will be merged.
#' @param unsorted.forward.counts Character vector of unsorted (input) count file paths generated by vl_ORFtag_pipeline().
#' Counts from all provided files will be merged.
#' @param unsorted.reverse.counts Character vector of insorted (input), REVERSED count file paths generated by vl_ORFtag_pipeline().
#' Counts from all provided files will be merged.
#' @param genome Genome annotation to be used. For now, only "mm10" and "hg38" are supported, and the corresponding non-first-exon .gtf will be used for sorting.
#' See at the bottom of the ?vl_ORFtag_pipeline help page to see how to generate custom non-first-exon start .gtf files for other genomes.
#' @param name Name prefix appended to output file names.
#' @param padj.cutoff The p.adjust cutoff to be used to call hits (>=). Default= 0.05.
#' @param log2OR.cutoff The log2OR cutoff to be used to call hits (>=). Default= 1.
#' @param log2OR.pseudocount The pseudocount to be added to avoid infinite values.
#' Note that only the log2OR will be affected, not the fisher p.values. Default= 1
#' @param binom.pseudocount The pseudocount to be added before binomial test (does not tolerate 0s). Default= 1
#' @param binom.padj.cutoff The p.adjust cutoff to be used to call binom_hits (<=). Default= 0.001.
#' @param output.suffix Name suffix appended to output file names. Default= "_vs_unsorted.txt".
#' @param output.folder.FC.file Output folder for FC files. Default= "db/FC_tables/ORFtag".
#'
#' @return Returns FC tables containing DESeq2-like columns.
#'
#' @examples
#' library(vlfunctions)
#' 
#' # 63 hits should be called with this dataset ----
#' vl_ORFtrap_call_hits_strandBias(sorted.forward.counts = c("db/exon_assignment/ORFtag/Activator2_sort_rep1_same_strand.txt",
#'                                                           "db/exon_assignment/ORFtag/Activator2_sort_rep2_same_strand.txt"),
#'                                 sorted.reverse.counts = c("db/exon_assignment/ORFtag/Activator2_sort_rep1_rev_strand.txt",
#'                                                           "db/exon_assignment/ORFtag/Activator2_sort_rep2_rev_strand.txt"),
#'                                 unsorted.forward.counts = c("db/exon_assignment/ORFtag/Activator2_input_rep1_same_strand.txt",
#'                                                             "db/exon_assignment/ORFtag/Activator2_input_rep2_same_strand.txt"),
#'                                 unsorted.reverse.counts = c("db/exon_assignment/ORFtag/Activator2_input_rep1_rev_strand.txt",
#'                                                             "db/exon_assignment/ORFtag/Activator2_input_rep2_rev_strand.txt"),
#'                                 genome = "mm10",
#'                                 name = "Activator_2",
#'                                 output.suffix = "_vs_input_strandBias.txt")
#'                                 
#' @export
vl_ORFtrap_call_hits_strandBias <- function(sorted.forward.counts, 
                                            sorted.reverse.counts, 
                                            unsorted.forward.counts,
                                            unsorted.reverse.counts,
                                            genome,
                                            name,
                                            padj.cutoff= 0.05,
                                            log2OR.cutoff= 1,
                                            log2OR.pseudocount= 1,
                                            binom.pseudocount= 1,
                                            binom.padj.cutoff= 0.001,
                                            output.suffix= "_vs_revStrand",
                                            output.folder.FC.file= "db/FC_tables/ORFtag")
{
  require(rtracklayer)
  require(data.table)
  require(GenomicRanges)
  
  # Checks
  if(!is.character(name) | length(name)!=1)
    stop("name should be a unique character specifying the name of the screen")
  if(anyDuplicated(sorted.forward.counts))
    stop("Duplicated filenames in sorted.forward.counts")
  if(anyDuplicated(unsorted.forward.counts))
    stop("Duplicated filenames in unsorted.forward.counts")
  if(!is.character(output.suffix) | length(output.suffix)!=1)
    stop("output.suffix should be a unique character specifying the name of the screen")
  
  # Import counts
  inputFw <- rbindlist(lapply(unsorted.forward.counts, fread))
  inputRev <- rbindlist(lapply(unsorted.reverse.counts, fread))
  sampleFw <- rbindlist(lapply(sorted.forward.counts, fread))
  sampleRev <- rbindlist(lapply(sorted.reverse.counts, fread))
  
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
  gtf <- switch(genome,
                "mm10" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
                "hg38" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_hg38.gtf")
  genes <- rtracklayer::import(gtf)
  genes <- as.data.table(mcols(genes)[, c("gene_id", "gene_name")])
  genes <- unique(genes)
  setorderv(genes, "gene_id")
  
  # Format dat
  dat <- merge(genes, dat, by= "gene_id", all.x= T)
  
  # Fisher test fw vs rev
  dat[count.sample.fw>=3, c("OR", "pval"):= {
    .t <- matrix(c(count.sample.fw,
                   count.sample.rev,
                   count.input.fw,
                   count.input.rev),
                 ncol= 2)
    .(fisher.test(.t, alternative = "greater")$estimate,
      fisher.test(.t+log2OR.pseudocount, alternative = "greater")$p.value)
  }, .(count.sample.fw, count.sample.rev, count.input.fw, count.input.rev)]
  dat[count.sample.fw>=3, log2OR:= log2(OR)]
  dat[count.sample.fw>=3, padj:= p.adjust(pval, "fdr")]
  dat[, hit:= padj<padj.cutoff & log2OR>=log2OR.cutoff]
  
  # Add binomial test
  dat[count.sample.fw>=3, binom_pval:= {
    binom.test(x = count.sample.fw+binom.pseudocount,
               n = count.sample.fw+count.sample.rev+binom.pseudocount,
               p = (count.input.fw+binom.pseudocount)/(count.input.fw+count.input.rev+binom.pseudocount),
               alternative = "greater")["p.value"]
  }, .(count.sample.fw, count.sample.rev, count.input.fw, count.input.rev)]
  dat[count.sample.fw>=3, binom_padj:= p.adjust(binom_pval, "fdr")]
  dat[, binom_hit:= binom_padj<binom.padj.cutoff]
  
  # Clean and save
  dat <- genes[dat, on= c("gene_id", "gene_name")]
  dat$OR <- dat$pval <- dat$binom_pval <- NULL
  dir.create(output.folder.FC.file, showWarnings = F, recursive = T)
  FC_table <- paste0(output.folder.FC.file, "/", name, output.suffix)
  fwrite(dat, 
         FC_table, 
         col.names = T, 
         row.names = F, 
         sep= "\t",
         quote= F, 
         na= NA)
  cat(paste0(name, ": ", sum(dat$hit, na.rm = T), " hits were called!\nFC file -> ", FC_table, "\n"))
}