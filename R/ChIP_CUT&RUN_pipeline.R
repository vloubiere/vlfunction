#' CUTNRUN pipeline
#' 
#' To use this this pipeline, please install the vl_function package using install_github("vloubiere/vlfunction")
#' in your interactive Rstudio session AND in the local installation of R 4.3.0 ('/software/f2022/software/r/4.3.0-foss-2022b/bin/R') 
#' which will be used to submit R sub-script.
#' 
#' The pipeline is split into two main functions. The first vl_CUTNRUN_processing() function aligns the reads and filters confident alignments.
#' It takes as input a (correctly formatted) metadata file, saves the processed metadata file and returns and/or submit the command lines to: \cr
#' 1/ extract reads from VBC .bam file. Output .fq files are saved tmp_folder/fq/\cr
#' 2/ trim the reads. Output .fq files are saved in tmp_folder/fq/ \cr
#' 3/ aligns them to mouse/human genome (see 'genome' column of the metadata table). Output .bam files are saved in tmp_folder/bam/ \cr
#' 4/ return alignment statistics. Output .txt files are saved in alignment_stats_output_folder/ \cr
#' 5/ return mapq30_statistics (confidently aligned reads). Output .txt files are saved in alignment_stats_output_folder/ \cr
#' 
#' The second function, vl_CUTNRUN_peakCalling(), takes care of the peak calling, and returns and/or submit the command lines to: \cr
#' 1/ peaks for each replicate. Output .narrowpeak files are saved in peaks_output_folder/ \cr
#' 2/ .bw tracks for each replicate. Output .bw files are saved in bw_output_folder/ \cr
#' 3/ Peaks called on the merged replicates. Output .narrowpeak files are saved in peaks_output_folder/ \cr
#' 4/ .bw tracks using merged replicates. Output .bw files are saved in bw_output_folder/ \cr
#' 5/ Confident peaks that are detected using merged reads but also each individual replicate. Output .narrowpeak files are saved in peaks_output_folder/ \cr
#'
#' @param metadata The path to a correctly formatted .xlsx metadata file or a data.table. See the template at '/groups/stark/vloubiere/projects/vl_pipelines/Rdata/metadata_CutNRun.xlsx'.
#' @param processed_metadata_output An .rds path where to save the processed metadata file, which contains the paths of all output files and will be used to manage them.
#' By default, when importing the metadata from an excel sheet, "_processed.rds" will be appended to the excel file path. 
#' @param alignment_stats_output_folder Output folder for alignment statistics. Default= "db/alignment_stats/CUTNRUN/".
#' @param peaks_output_folder Output folder for peak files. Default= "db/peaks/CUTNRUN/".
#' @param bw_output_folder  Output folder for .bw tracks. Default= "db/bw/CUTNRUN/".
#' @param tmp_folder Output folder for temporary files (.fq, .bam). Default= "/scratch/stark/vloubiere/ORFtag".
#' @param cores Number of cores per job. Default= 8.
#' @param mem Memory per job (in Go). Default= 32.
#' @param overwrite Should existing files be overwritten?
#' @param submit Should the command be submitted? Default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs/CUTNRUN/processing".
#' @param time The time required for the SLURM scheduler. Default= '1-00:00:00'.
#'
#' @return Return a data.table containing, for each sampleID, the concatenated commands that are required to process the data.
#' These commands can then be submitted either directly via the function, or using vl_bsub()....
#'
#' @examples
#' # Process example dataset ----
#' library(vlfunctions)
#' meta <- vl_metadata_CUTNRUN
#' vl_CUTNRUN_processing(metadata = meta,
#'                       processed_metadata_output = "Rdata/metadata_CutNRun_processed.rds",
#'                       cores= 8,
#'                       mem= 64,
#'                       overwrite= FALSE,
#'                       submit= TRUE)
#'                       
#' # Once the jobs have finished to work, you can check the size of output files:
#' processed <- readRDS("Rdata/metadata_CutNRun_processed.rds")
#' file.size(unlist(unique(processed[, fq1:mapq30_stats])))
#' 
#' # Peak calling (see ?vl_CUTNRUN_peakCalling)----
#' vl_CUTNRUN_peakCalling(processed_metadata = processed, # You can also provide the path to the file.
#'                        extsize = 300,
#'                        cores = 8,
#'                        mem = 64,
#'                        overwrite = FALSE,
#'                        submit = TRUE)
#'
#' @export
vl_CUTNRUN_processing <- function(metadata, ...) UseMethod("vl_CUTNRUN_processing")

#' @describeIn vl_CUTNRUN_processing for excel files path
#' @export
vl_CUTNRUN_processing.character <- function(metadata,
                                            processed_metadata_output= gsub(".xlsx$", "_processed.rds", metadata),
                                            ...)
{
  sheet <- readxl::read_xlsx(metadata)
  start <- cumsum(sheet[[1]]=="user")>0
  sheet <- sheet[(start),]
  meta <- as.data.table(sheet[-1,])
  names(meta) <- unlist(sheet[1,])
  vl_CUTNRUN_processing(metadata= meta,
                        processed_metadata_output= processed_metadata_output,
                        ...)
}

#' @describeIn vl_CUTNRUN_processing default method
#' @export
vl_CUTNRUN_processing.default <- function(metadata,
                                          processed_metadata_output,
                                          alignment_stats_output_folder= "db/alignment_stats/CUTNRUN/",
                                          peaks_output_folder= "db/peaks/CUTNRUN/",
                                          bw_output_folder= "db/bw/CUTNRUN/",
                                          tmp_folder= "/scratch/stark/vloubiere",
                                          cores= 8,
                                          mem= 64,
                                          overwrite= FALSE,
                                          submit= FALSE,
                                          wdir= getwd(),
                                          logs= "db/logs/CUTNRUN/processing",
                                          time= '1-00:00:00')
{
  # Import metadata and check format ----
  meta <- data.table::copy(metadata)
  cols <- c("user","method","experiment","antibody","cell_line","treatment","replicate","condition","sampleID","barcode","i5","genome","input","layout","sequencer","bam_path")
  if(!all(cols %in% names(meta)))
    stop(paste("Columns missing ->", paste(cols[!cols %in% names(meta)], collapse = "; ")))
  if(any(!na.omit(meta$input) %in% meta$sampleID))
    stop("Some input(s) have no correspondence in the sampleID(s) column!")
  if(!all(meta[, condition==paste0(antibody, "_", cell_line, "_", treatment, "_", experiment)]))
    warning("To avoid confusion between samples, it would be safer if condition was the concatenation of ab, cell_line, treatment and experiment")
  if(!all(meta[, sampleID==paste0(antibody, "_", cell_line, "_", treatment, "_", experiment, "_", replicate)]))
    warning("sampleID should be the concatenation of ab, cell_line, treatment, experiment and replicate (separated by '_'). Any samples with the same sampleID will be collapsed!")
  if(!all(meta$layout %in% c("SINGLE", "PAIRED")))
    stop("layout column should only contain 'SINGLE' or 'PAIRED'")
  
  # Check bowtie2 idx ----
  if(any(!meta$genome %in% c("mm10", "hg38")))
    stop("Only mm10 and hg38 are supported! For other genomes, please provide path to the corresponding bowtie 2 index.")
  
  # Check fq files ----
  if("fq1" %in% names(meta) && any(grepl(".fq.gz$", na.omit(meta$fq1))))
    stop("If provided, fq1 files should be gzipped and their names should end with '.fq.gz'")
  if("fq2" %in% names(meta) && any(grepl(".fq.gz$", na.omit(meta$fq2))))
    stop("If provided, fq2 files should be gzipped and their names should end with '.fq.gz'")
  
  # Generate output paths ----
  meta[is.na(fq1) & !is.na(bam_path), fq1:=
         paste0(tmp_folder, "/CUTNRUN/fq/", gsub(".bam", "", basename(bam_path)), "_", make.unique(sampleID),
                fifelse(layout=="PAIRED", "_1.fq.gz", ".fq.gz"))]
  meta[is.na(fq2) & layout=="PAIRED", fq2:= gsub("_1.fq.gz$", "_2.fq.gz", fq1)]
  meta[, fq1_trim:= gsub(".fq.gz$", fifelse(layout=="PAIRED", "_val_1.fq.gz", "_trimmed.fq.gz"), fq1), layout]
  meta[, fq2_trim:= gsub(".fq.gz$", "_val_2.fq.gz", fq2)]
  # re-sequencing are merged from this step on!
  meta[, bam:= paste0(tmp_folder, "/CUTNRUN/bam/", sampleID, ".bam")]
  # Join to retrieve input bam
  meta[meta, bam_input:= i.bam, on= "input==sampleID"] # Input bam
  meta[, alignment_stats:= paste0(alignment_stats_output_folder, gsub(".bam$", "_stats.txt", basename(bam)))]
  meta[, mapq30_stats:= gsub("_stats.txt$", "_mapq30_stats.txt", alignment_stats)]
  meta[, peaks_file:= paste0(peaks_output_folder, "replicates/", sampleID, "_peaks.narrowPeak")]
  meta[, bw_file:= paste0(bw_output_folder, "replicates/", sampleID, ".bw")]
  # Grouping per condition to check N replicates
  meta[, twoReps:= length(unique(bam))>1, condition] 
  meta[(twoReps), merged_peaks_file:= paste0(peaks_output_folder, "merge/", condition, "_merged_peaks.narrowPeak")]
  meta[(twoReps), merged_bw_file:= paste0(bw_output_folder, "merge/", condition, "_merged.bw")]
  meta[(twoReps), confident_peaks_file:= paste0(peaks_output_folder, "confident/", condition, "_confident_peaks.narrowPeak")]
  
  # Save processed metadata ----
  if(!grepl(".rds$", processed_metadata_output))
    stop("processed_metadata_output should end up with a .rds extension") else
      saveRDS(meta, processed_metadata_output)
  
  # Create output directories ----
  tmpSort <- paste0(tmp_folder, "/ORFtag/bam/sort/") # tmp sort bam files
  dirs <- c(logs,
            tmpSort,
            na.omit(unique(dirname(unlist(meta[, fq1:confident_peaks_file])))))
  if(any(!dir.exists(dirs)))
  {
    sapply(dirs, dir.create, showWarnings = F, recursive = T)
    print("Output directories were created in db/ and scratch/!")
  }
  
  # Print conditions ----
  if(any(meta$twoReps))
  {
    cat(paste(length(unique(meta$condition)), "conditions detected, of which", length(unique(meta[(twoReps), condition])), "had >1 replicates:\n"))
    meta[(twoReps), cat(paste0("> ", condition, " -> ", paste0(unique(sampleID), collapse= "; "), "\n")), condition]
  }
  if(any(!meta$twoReps))
  {
    cat(paste(length(unique(meta[(!twoReps), condition])), "samples had only one replicate:\n"))
    meta[(!twoReps), cat(paste0("> ", condition, " -> ", paste0(unique(sampleID), collapse= "; "), "\n")), condition]
  }
  
  # Demultiplex VBC bam file ----
  meta[, demultiplex_cmd:= {
    if(.N!=1)
      stop("Unique bam file should be provided!")
    if(overwrite | !file.exists(fq1))
    {
      demultScript <- system.file("CUTNRUN_pipeline",
                                  ifelse(layout=="PAIRED", "demultiplex_pe.pl", "demultiplex_se.pl"),
                                  package = "vlfunctions")
      fq_prefix <- gsub("_1.fq.gz$|.fq.gz$", "", fq1)
      cmd <- paste("samtools view -@", cores-1, bam_path,
                   "| perl ", demultScript,
                   # "| head -n 1000000 | perl ", demultScript, # For tests
                   paste0("^B2:Z:", i5),
                   paste0("^BC:Z:", barcode),
                   fq_prefix,
                   "; gzip -f", gsub(".gz$", "", fq1))
      if(!is.na(fq2))
        cmd <- paste0(cmd, "; gzip -f ", gsub(".gz$", "", fq2))
      cmd
    }
  }, .(bam_path, layout, i5, barcode, fq1, fq2)]
  
  # Trim illumina adapters ----
  meta[, trim_cmd:= {
    if(overwrite | !file.exists(fq1_trim) | (!is.na(fq2_trim)  && !file.exists(fq2_trim)))
    {
      if(layout=="SINGLE")
        paste0("trim_galore --gzip -o ", dirname(fq1_trim), "/ ", fq1) else if(layout=="PAIRED")
          paste0("trim_galore --gzip --paired -o ", dirname(fq1_trim), "/ ", fq1, " ", fq2)
    }
  }, .(layout, fq1, fq1_trim, fq2, fq2_trim)]
  
  # Alignment ----
  meta[, aln_cmd:= {
    if(overwrite | !file.exists(bam))
    {
      # BOWTIE 2 ----
      inputFiles <- if(layout=="SINGLE")
        paste("-U", paste0(unique(fq1_trim), collapse= ",")) else if(layout=="PAIRED")
          paste("-1", paste0(unique(fq1_trim), collapse= ","),
                "-2", paste0(unique(fq2_trim), collapse= ","))
      x <- switch(genome,
                  "mm10" = "/groups/stark/vloubiere/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome",
                  "hg38" = "/groups/stark/vloubiere/genomes/Homo_sapiens/hg38/Bowtie2Index/genome")
      paste("bowtie2 -p", cores,
            "-x", x,
            inputFiles,
            "2>", alignment_stats, # Return alignment statistics
            "| samtools sort -@", cores-1, "-T", tmpSort, # Sort with temp files in scratch folder
            "- | samtools view -@", cores-1, "-b -q 30 -o", bam, # Filter
            "; samtools stats", bam, "-@", cores-1, "| grep ^SN>", mapq30_stats) # Add filtering statistics
    }
  }, .(layout, bam, genome, alignment_stats, mapq30_stats)]

  # Return commands ----
  load_cmd <- paste(c(paste("cd", wdir),
                      "module load build-env/2020",
                      "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                      "module load samtools/1.9-foss-2018b",
                      "module load bowtie2/2.3.4.2-foss-2018b"), collapse = "; ")
  cols <- intersect(c("demultiplex_cmd", "trim_cmd", "aln_cmd"),
                    names(meta))
  cmd <- if(length(cols))
    meta[, .(cmd= paste0(c(load_cmd, unique(na.omit(unlist(.SD)))),
                         collapse = "; ")), sampleID, .SDcols= cols] else
      data.table()
  
  # Submit commands ----
  if(nrow(cmd))
  {
    if(submit)
      cmd[, {
        vl_bsub(cmd, 
                cores= cores, 
                m = mem, 
                name = "CUTNRUN", 
                t = time,
                o= paste0(normalizePath(logs), "/", sampleID),
                e= paste0(normalizePath(logs), "/", sampleID))
      }, .(sampleID, cmd)] else
        return(cmd)
  }else
    warning("All output files already existed! No command submitted ;). Consider overwrite= T if convenient.")
}

#' Title
#'
#' @param processed_metadata Path to the metadata file generated by vl_CUTNRUN_processing, or the corresponding data.table.
#' @param Rpath The path to the Rscript executable to use, on which the latest version of the vlfunction package should be installed. Default= /software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript
#' @param extsise The extsize to be used, meaning that the building of the peak model will be skipped (--nomodel). To build the model, set extsize= NA. Default= 300. 
#' @param cores Number of cores per job. Default= 8.
#' @param mem Memory per job (in Go). Default= 32.
#' @param overwrite Should existing files be overwritten?
#' @param submit Should the command be submitted? default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs/CUTNRUN/peak_calling".
#' @param time The time required for the SLURM scheduler. Default= '1-00:00:00'.
#'
#' @return Command lines for peak calling.
#'
#' @examples
#' # Peak calling example dataset ----
#' library(vlfunctions)
#' 
#' # Check all processed files exist
#' processed <- readRDS("Rdata/metadata_CutNRun_processed.rds")
#' file.exists(na.omit(unlist(processed[, fq1:mapq30_stats])))
#' 
#' vl_CUTNRUN_peakCalling(processed_metadata = processed, # You could also just provide the path, "Rdata/metadata_CutNRun_processed.rds".
#'                        extsize = 300,
#'                        cores = 8,
#'                        mem = 64,
#'                        overwrite = FALSE,
#'                        submit = TRUE)
#' @export
vl_CUTNRUN_peakCalling <- function(processed_metadata, ...) UseMethod("vl_CUTNRUN_peakCalling")

#' @describeIn vl_CUTNRUN_peakCalling for processed_metadata file path
#' @export
vl_CUTNRUN_peakCalling.character <- function(processed_metadata, ...)
{
  meta <- readRDS(processed_metadata)
  vl_CUTNRUN_peakCalling(processed_metadata= meta, ...)
}

#' @describeIn vl_CUTNRUN_processing default method
#' @export
vl_CUTNRUN_peakCalling.default <- function(processed_metadata,
                                           Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                           extsize= 300,
                                           cores= 8,
                                           mem= 64,
                                           overwrite= FALSE,
                                           submit= FALSE,
                                           wdir= getwd(),
                                           logs= "db/logs/CUTNRUN/peak_calling",
                                           time= '1-00:00:00')
{
  # Hard copy metadata ----
  meta <- data.table::copy(processed_metadata)
  
  # Check bowtie2 idx ----
  if(any(!meta$genome %in% c("mm10", "hg38")))
    stop("Only mm10 and hg38 are supported! For other genomes, please provide path to the corresponding bowtie 2 index.")
  
  # Check if at least two replicates and input specified ----
  if(any(is.na(meta$input)))
    meta[is.na(input), print(paste("No input found for", sampleID, "-> peak calling without input")), sampleID]
  if(any(!meta$twoReps))
    meta[(!twoReps), print(paste("No replicates found for", sampleID, "-> confident peaks calling skipped")), sampleID]
  
  # Create missing output dirs ----
  if(!dir.exists(logs))
  {
    dir.create(logs, showWarnings = F, recursive = T)
    print(paste(logs, "directory was created!"))
  }
  
  # Peak calling on individual replicates ----
  meta[, peak_cmd:= {
    if(overwrite | !file.exists(peaks_file))
    {
      # Genome
      gAbb <- fcase(genome=="mm10", "mm",
                    genome=="hg38", "hg")
      # Check input
      .input <- ifelse(is.na(bam_input), "", paste("-c", bam_input))
      # Command
      paste("macs2 callpeak -B --SPMR --keep-dup 1 -g", gAbb, 
            ifelse(layout=="PAIRED", "-f BAMPE", ""),
            ifelse(is.na(extsize), "", paste("--nomodel --extsize", extsize)),
            "--outdir", paste0(dirname(peaks_file), "/"),
            "-t", bam,
            .input,
            "-n", gsub("_peaks.narrowPeak", "", basename(peaks_file)))
            # if(broad) # Broad signal (K27me3...)
            #   cmd <- paste0(cmd, " --broad")
    }
  }, .(layout, genome, bam, bam_input, peaks_file)]
  
  # Bedgraph to bigwig individual rep ----
  meta[, bw_cmd:= {
    if(overwrite | !file.exists(bw_file))
    {
      # Usage: Rscript script.R input_file.bdg output_file.bw genome
      paste(Rpath, system.file("CUTNRUN_pipeline", "bedgraph_to_bigwig.R", package = "vlfunctions"),
            gsub("_peaks.narrowPeak", "_treat_pileup.bdg", peaks_file),
            bw_file,
            genome)
    }
  }, .(genome, bam, bam_input, peaks_file, bw_file)]
  
  # Merged peaks ----
  meta[(twoReps), merged_peak_cmd:= {
    if(overwrite | !file.exists(merged_peaks_file))
    {
      # Genome
      gAbb <- fcase(genome=="mm10", "mm",
                    genome=="hg38", "hg")
      # Check input
      .input <- paste(unique(na.omit(bam_input)), collapse= " ")
      .input <- ifelse(length(.input), paste("-c", .input), "")
      # Command
      paste("macs2 callpeak -B --SPMR --keep-dup 1 -g", gAbb,
            ifelse(layout=="PAIRED", "-f BAMPE", ""),
            ifelse(is.na(extsize), "", paste("--nomodel --extsize", extsize)),
            "--outdir", paste0(dirname(merged_peaks_file), "/"),
            "-t", paste0(unique(bam), collapse= " "),
            .input,
            "-n", gsub("_peaks.narrowPeak", "", basename(merged_peaks_file)))
    }
  }, .(genome, condition, merged_peaks_file)]
  
  # Merged bigwig ----
  meta[(twoReps), merged_bw_cmd:= {
    if(overwrite | !file.exists(merged_bw_file))
    {
      # Usage: Rscript script.R input_file.bdg output_file.bw genome
      paste(Rpath, system.file("CUTNRUN_pipeline", "bedgraph_to_bigwig.R", package = "vlfunctions"),
            gsub("_peaks.narrowPeak", "_treat_pileup.bdg", merged_peaks_file),
            merged_bw_file,
            genome)
    }
  }, .(genome, condition, merged_peaks_file, merged_bw_file)]
  
  # Confident peaks ----
  meta[(twoReps), conf_peak_cmd:= {
    if(overwrite | !file.exists(confident_peaks_file))
    {
      paste(Rpath, system.file("CUTNRUN_pipeline", "confident_peaks.R", package = "vlfunctions"),
            paste0(unique(peaks_file), collapse= ","),
            merged_peaks_file,
            confident_peaks_file)
    }
  }, .(condition, merged_peaks_file, confident_peaks_file)]
  
  # Return commands ----
  load_cmd <- paste(c(paste("cd", wdir),
                      "module load build-env/2020",
                      "module load macs2/2.1.2.1-foss-2018b-python-2.7.15"),
                    collapse = "; ")
  cols <- intersect(c("peak_cmd", "bw_cmd", "merged_peak_cmd", "merged_bw_cmd", "conf_peak_cmd"),
                    names(meta))
  cmd <- if(length(cols))
    meta[, .(cmd= paste0(c(load_cmd, unique(na.omit(unlist(.SD)))),
                         collapse = "; ")), condition, .SDcols= cols] else
                           data.table()
  
  # Submit commands ----
  if(nrow(cmd))
  {
    if(submit)
      cmd[, {
        vl_bsub(cmd, 
                cores= cores, 
                m = mem, 
                name = "CUTNRUN", 
                t = time,
                o= paste0(normalizePath(logs), "/", condition),
                e= paste0(normalizePath(logs), "/", condition))
      }, .(condition, cmd)] else
        return(cmd)
  }else
    warning("All output files already existed! No command submitted ;). Consider overwrite= T if convenient.")
}