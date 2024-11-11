#' RNA-seq pipeline
#' 
#' To use this this pipeline, please install the vl_function package using install_github("vloubiere/vlfunction")
#' in your interactive Rstudio session AND in the local installation of R 4.3.0 ('/software/f2022/software/r/4.3.0-foss-2022b/bin/R') 
#' which will be used to submit R sub-scripts.
#' 
#' The pipeline is split into two main function The first function vl_RNAseq_processing() demultiplexes tar files, aligns the reads and filters confident alignments.
#' It takes as input a (correctly formatted) metadata file, saves the processed metadata file and returns and/or submit the command lines to: \cr
#' 1/ extract reads from VBC .tar.gz file. Output .fq.gz files are saved in tmp_folder/fq/ \cr
#' 2/ trim the reads. Output .fq.gz files are saved in tmp_folder/fq/ \cr
#' 3/ align the reads to the provided genome (see 'genome' column of the metadata table, only mm10 is supported for now). Output .bam files are saved in tmp_folder/bam/ and alignment statistics in the alignment_stats_output_folder/ \cr
#' 4/ compute read counts per exon, and returns the total number of reads per gene_id. Counts statistics will be save in alignment_stats_output_folder/ (together with alignment stats) and count_tables in count_table_output_folder/  \cr
#' 5/ output .bw tracks that are stored in bw_output_folder/
#' 
#' The second function, vl_RNAseq_DESeq2() can be used for differential analysis. It returns: \cr
#' 1/ a .dds DESeq2 object for each experiment. Output .dds (rds) files are stored in the dds_output_folder/ \cr
#' 2/ the fold-change tables corresponding to the comparison of 'DESeq2_condition' versus 'DESeq2_control' conditions (see colnames of the processed_metadata file). Output .txt files are stored in the FC_output_folder/ \cr
#' 3/ A barplot showing the statistics of each dds file, which will be saved in the PDF_output_folder. \cr
#' 4/ MA plots corresponding to each comparison, which will be saved in the PDF_output_folder. \cr
#'
#' @param metadata The path to a correctly formatted .xlsx metadata file or a data.table. See the template at '/groups/stark/vloubiere/projects/vl_pipelines/Rdata/metadata_RNAseq.xlsx'.
#' @param processed_metadata_output An .rds path where to save the processed metadata file, which contains the paths of all output files and will be used to manage them.
#' By default, when importing the metadata from an excel sheet, "_processed.rds" will be appended to the excel file path. 
#' @param count_table_output_folder Output folder for count table files. Default= "db/count_tables/RNAseq/".
#' @param alignment_stats_output_folder Output folder for alignment statistics. Default= "db/alignment_stats/RNAseq/".
#' @param bw_output_folder Output folder for .bw tracks. Default= "db/bw/RNAseq/".
#' @param tmp_folder Output folder for temporary files (.fq, .bam). Default= "/scratch/stark/vloubiere".
#' @param Rpath Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of cores per job. Default= 8.
#' @param mem Memory per job (in Go). Default= 64.
#' @param overwrite Should existing files be overwritten?
#' @param counts.exists If the output count_table file already exists and overwrite is set to FALSE, skips demultiplexing and alignment. Default= TRUE.
#' @param submit Should the command be submitted? Default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs/CUTNRUN/processing".
#' @param time The time required for the SLURM scheduler. Default= '2-00:00:00'.
#'
#' @return Return a data.table containing, for each sampleID, the concatenated commands that are required to process the data.
#' These commands can then be submitted either directly via the function, or using vl_bsub()...
#'
#' @examples
#' # Before starting:
#' # Update the local version of the vlfunctions package
#' devtools::install_github("vloubiere/vlfunction")
#' 
#' # Chose which local R executable to use
#' localRPath <- "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"
#' # Make sure that the last version of the vlfunctions package is installed
#' system(paste(localRPath, system.file("update_package_from_git.R", package = "vlfunctions")))
#' 
#' # Process example dataset ----
#' example <- vl_metadata_PROseq
#' vl_PROseq_processing(metadata = example,
#'                      processed_metadata_output = "Rdata/metadata_PROseq_processed.rds",
#'                      Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript", 
#'                      cores= 8,
#'                      mem= 64,
#'                      overwrite= FALSE,
#'                      submit= FALSE)
#' 
#' # Once the jobs have finished to work, you can check the size of output files:
#' processed <- readRDS("Rdata/metadata_PROseq_processed.rds")
#' file.size(na.omit(unlist(processed[, fq1:bw])))
#' 
#' # Differential analysis
#' See ?vl_RNAseq_DESeq2()
#'
#' @export
vl_RNAseq_processing <- function(metadata, ...) UseMethod("vl_RNAseq_processing")

#' @describeIn vl_RNAseq_processing for excel files path
#' @export
vl_RNAseq_processing.character <- function(metadata,
                                           processed_metadata_output= gsub(".xlsx$", "_processed.rds", metadata),
                                           ...)
{
  sheet <- readxl::read_xlsx(metadata)
  start <- cumsum(sheet[[1]]=="user")>0
  sheet <- sheet[(start),]
  meta <- as.data.table(sheet[-1,])
  names(meta) <- unlist(sheet[1,])
  vl_RNAseq_processing(metadata= meta,
                       processed_metadata_output= processed_metadata_output,
                       ...)
}

#' @describeIn vl_RNAseq_processing default method
#' @export
vl_RNAseq_processing.default <- function(metadata,
                                         processed_metadata_output,
                                         alignment_stats_output_folder= "db/alignment_stats/RNAseq/",
                                         count_table_output_folder= "db/count_tables/RNAseq/",
                                         bw_output_folder= "db/bw/RNAseq/",
                                         tmp_folder= "/scratch/stark/vloubiere",
                                         Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                         cores= 8,
                                         mem= 64,
                                         overwrite= FALSE,
                                         counts.exists= TRUE,
                                         submit= FALSE,
                                         wdir= getwd(),
                                         logs= "db/logs/RNAseq/processing",
                                         time= '2-00:00:00')
{
  # Import metadata and check format ----
  meta <- data.table::copy(metadata)
  cols <- c("user", "batch", "experiment", "cell_line", "treatment", "replicate", "sampleID",
            "DESeq2_name", "DESeq2_condition", "DESeq2_control", "barcode", "i5", "genome",
            "layout", "sequencer", "fastq_path")
  if(!all(cols %in% names(meta)))
    stop(paste("Columns missing ->", paste(cols[!cols %in% names(meta)], collapse = "; ")))
  if(any(!na.omit(meta$DESeq2_control) %in% meta$DESeq2_condition))
    stop("Some DESeq2_control values have no correspondence in the DESeq2_condition column!")
  if(!all(meta[, sampleID==paste0(cell_line, "_", treatment, "_", experiment, "_", replicate)]))
    warning("sampleID should be the concatenation of cell_line, treatment, experiment and replicate (separated by '_'). Any samples with the same sampleID will be collapsed!")
  if(!all(meta$layout %in% c("SINGLE", "PAIRED")))
    stop("layout column should only contain 'SINGLE' or 'PAIRED'")
  if(any(meta$layout!="PAIRED"))
    stop("Single-end reads are not implemented yet :(!")
  if(any(meta[, length(unique(layout)), .(sampleID, genome)]$V1>1))
    stop("For unique sampleID / genome combination(s), the layout (paired-end, single-end) should always be consistent!")
  
  # Check bowtie2 idx ----
  if(any(!meta$genome %in% c("mm10", "hg38")))
    stop("Only mm10 and hg38 are supported! For other genomes, please provide path to the corresponding bowtie 2 index.")
  if(any(!meta$spikein_genome %in% c("dm3")))
    stop("Only dm3 is supported for spike-in! For other genomes, please provide path to the corresponding bowtie 2 index.")
  
  # Generate output paths ----
  meta[, fq1:= paste0(tmp_folder, "/RNAseq/fq/", gsub(".tar.gz", paste0("_", barcode, "_", i5, "_1.fq.gz"), basename(fastq_path))), .(fastq_path, barcode, i5)]
  meta[layout=="PAIRED", fq2:= gsub("_1.fq.gz$", "_2.fq.gz", fq1)]
  meta[, fq1_trim:= gsub(".fq.gz$", fifelse(layout=="PAIRED", "_val_1.fq.gz", "_trimmed.fq.gz"), fq1), layout]
  meta[, fq2_trim:= gsub(".fq.gz$", "_val_2.fq.gz", fq2)]
  # re-sequencing are merged from this step on!
  meta[, bam:= paste0(tmp_folder, "/RNAseq/bam/", sampleID, "_", genome, ".bam")]
  # reads statistics
  meta[, read_stats:= paste0(alignment_stats_output_folder, experiment, "/", sampleID, "_", genome, "_statistics.txt")]
  # Read counts
  meta[, count_table:= paste0(count_table_output_folder, experiment, "/", sampleID, "_", genome, "_counts.txt")]
  # # bw tracks
  meta[, bw:= paste0(bw_output_folder, experiment, "/", sampleID, ".bw")]
  
  # Save processed metadata ----
  if(!grepl(".rds$", processed_metadata_output))
    stop("processed_metadata_output should end up with a .rds extension") else
      saveRDS(meta, processed_metadata_output)
  
  # Create output directories ----
  dirs <- c(logs,
            na.omit(unique(dirname(unlist(meta[, fq1:bw])))))
  if(any(!dir.exists(dirs)))
  {
    sapply(dirs, dir.create, showWarnings = F, recursive = T)
    outDirs <- paste0(c(count_table_output_folder, bw_output_folder, tmp_folder), collapse= ", ")
    print(paste("Output directories were created in:", outDirs, "!"))
  }
  
  # Check whether demultiplexing and alignment should be performed ----
  meta[, preProcess:= ifelse(counts.exists, !file.exists(count_table), TRUE)] # Check if demultiplex/alignment should be performed
  
  # Extract fastq files ----
  meta[, extract_cmd:= {
    if(overwrite | (preProcess & any(!file.exists(c(fq1, fq2)))))
    {
      fq_prefix <- normalizePath(gsub("_1.fq.gz$", "", fq1), mustWork = F)
      paste("perl",
            fcase(layout=="PAIRED", system.file("RNAseq_pipeline", "RNAseq_demultiplexing_tar_pe.pl", package = "vlfunctions"),
                  layout=="SINGLE", system.file("RNAseq_pipeline", "RNAseq_demultiplexing_tar_se.pl", package = "vlfunctions")),
            fastq_path,
            barcode,
            i5,
            fq_prefix)
    }
  }, .(layout, fq1, fastq_path, barcode, i5, preProcess)]
  
  # Trim illumina adapters ----
  meta[, trim_cmd:= {
    if(overwrite | (preProcess & any(!file.exists(c(fq1_trim, fq2_trim)))))
    {
      if(layout=="SINGLE")
        paste0("trim_galore --gzip -o ", dirname(fq1_trim), "/ ", fq1) else if(layout=="PAIRED")
          paste0("trim_galore --gzip --paired -o ", dirname(fq1_trim), "/ ", fq1, " ", fq2)
    }
  }, .(layout, fq1, fq1_trim, fq2, fq2_trim, preProcess)]
  
  # UMI counts ----
  meta[, align_count_cmd:= {
    # [required] 1/ A comma-separated list of fq1 files \n
    # [required] 2/ A comma-separated list of fq2 files \n
    # [required] 3/ Is the data paired-end?
    # [required] 4/ subreadr index prefix \n
    # [required] 5/ Output bam path \n
    # [required] 6/ Path to the gtf file that was used to generate the subreadr index \n
    # [required] 7/ Output statistics file \n
    # [required] 8/ Output count table file \n"
    if(overwrite | !file.exists(count_table))
    {
      idx <- switch(genome,
                    "mm10"= "/groups/stark/vloubiere/genomes/Mus_musculus/subreadr_mm10/subreadr_mm10_index")
      gtf <- switch(genome,
                    "mm10"= "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/gencode.vM25.basic.annotation.gtf.gz")
      paste(Rpath,
            system.file("RNAseq_pipeline", "align_and_count.R", package = "vlfunctions"),
            paste0(fq1, collapse = ","),
            paste0(fq2, collapse = ","),
            layout=="PAIRED",
            idx,
            bam,
            gtf,
            read_stats,
            count_table)
    }
  }, .(genome, layout, bam, read_stats, count_table)]
  
  # Generate bw files ----
  meta[, bw_cmd:= {
    # Please specify:
    # [required] 1/ Input bam file
    # [required] 2/ Output bw file
    if(overwrite | !file.exists(bw))
    {
      paste(Rpath,
            system.file("RNAseq_pipeline", "generate_bw_files.R", package = "vlfunctions"),
            bam,
            bw)
    }
  }, .(bam, bw)]
  
  # Return commands ----
  load_cmd <- paste(c(paste("cd", wdir),
                      "module load build-env/2020",
                      "module load trim_galore/0.6.0-foss-2018b-python-2.7.15"),
                    collapse = "; ")
  cols <- intersect(c("extract_cmd", "trim_cmd", "align_count_cmd", "bw_cmd"),
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
                name = "RNAseq", 
                t = time,
                o= paste0(normalizePath(logs), "/", sampleID),
                e= paste0(normalizePath(logs), "/", sampleID))
      }, .(sampleID, cmd)] else
        return(cmd)
  }else
    warning("All output files already existed! No command submitted ;). Consider overwrite= T if convenient.")
}

#' DESeq2 analysis RNAseq pipeline
#'
#' @param processed_metadata Path to the metadata file generated by vl_PROseq_processing, or the corresponding data.table.
#' @param FC_metadata_output An .rds path where to save the metadata file, which contains the directories containing .dds (DESeq2 objects) and FC table files and will be used to manage them.
#' By default, when the processed_metadata is a path to a processed_metadata files, "_FC_tables.rds" will be appended to the excel file path. 
#' @param Rpath Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param count_table_output_folder Output folder for count tables. Default= "db/count_tables/PROseq/".
#' @param normalization A vector of normalization approaches to be used. Possible values are are "default" (DESeq2 default), "libSize", "spikeIn" and "combined" (Vanja's method: note that changing the control sample will change the outcome). Default= c("libSize", "spikeIn").
#' @param dds_output_folder Output folder for .dds files (DESeq2 objects). Default= "db/dds/PROseq/".
#' @param FC_output_folder Output folder for FC tables. Default= "db/FC_tables/PROseq/".
#' @param PDF_output_folder Output folder for .pdf files of MA plots and statistics. Default= "pdf/PROseq/".
#' @param cores Number of cores per job. Default= 4.
#' @param mem Memory per job (in Go). Default= 16.
#' @param overwrite Should the count tables be overwritten? Note that FC tables will always be overwritten. Default= FALSE
#' @param submit Should the command be submitted? Default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs/CUTNRUN/peak_calling".
#' @param time The time required for the SLURM scheduler. Default= '1-00:00:00'.
#'
#' @return Differential analyses using DESeq2.
#' @export
#'
#' @examples
#' # List annotation files for features of interest
#' annot <- data.table(genome= "mm10",
#'                     feature= c("promoter", "gene_body", "transcript"),
#'.                    annotation.file= c("/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_promoters.rds",
#'                                        "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_genebody.rds",
#'                                        "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_transcript.rds"))
#' 
#' # Import processed metadata file output by vl_PROseq_processing()
#' processed <- readRDS("Rdata/metadata_PROseq_processed.rds")
#' 
#' # Differential analysis
#' vl_PROseq_DESeq2(processed,
#'                  FC_metadata_output = "Rdata/metadata_PROseq_FC_tables.rds",
#'                  annotation.files= annot,
#'                  overwrite = FALSE,
#'                  submit = TRUE)
#' 
#' # Once the commands ran, you can use the FC metadata to retrieve the corresponding files:
#' FC <- readRDS("Rdata/metadata_PROseq_FC_tables.rds")
#' 
#' # For example, to retrieve spikeIn normalized files
#' FC[norm=="spikeIn"]
#' 
#' # Or for a specific feature
#' FC[feature=="promoter"]
#'                  
vl_RNAseq_DESeq2 <- function(processed_metadata, ...) UseMethod("vl_RNAseq_DESeq2")

#' @describeIn vl_RNAseq_DESeq2  for processed_metadata file path
#' @export
vl_RNAseq_DESeq2.character <- function(processed_metadata,
                                       FC_metadata_output= gsub("processed.rds$", "FC_tables.rds", processed_metadata),
                                       ...)
{
  meta <- readRDS(processed_metadata)
  vl_RNAseq_DESeq2(processed_metadata= meta,
                   FC_metadata_output= FC_metadata_output,
                   ...)
}

#' @describeIn vl_RNAseq_DESeq2 default method
#' @export
vl_RNAseq_DESeq2.default <- function(processed_metadata,
                                     FC_metadata_output,
                                     Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                     dds_output_folder= "db/dds/RNAseq/",
                                     FC_output_folder= "db/FC_tables/RNAseq/",
                                     PDF_output_folder= "pdf/RNAseq/",
                                     cores= 2,
                                     mem= 8,
                                     overwrite= FALSE,
                                     submit= FALSE,
                                     wdir= getwd(),
                                     logs= "db/logs/RNAseq/DESeq2",
                                     time= '01:00:00')
{
  # Hard copy metadata ----
  meta <- data.table::copy(processed_metadata)
  meta$bam_path <- meta$fastq_path <- meta$barcode <- meta$i5 <- meta$layout <- meta$sequencer <- NULL
  meta$fq1 <- meta$fq1_trim <- meta$fq2 <- meta$fq2_trim <- meta$bam <- NULL
  meta$bw <- NULL
  meta <- unique(meta) # Collapse re-sequencing
  
  # Check Number of replicates ----
  reps <- meta[, .N, .(experiment, DESeq2_condition)][N<2]
  if(nrow(reps))
  {
    reps[, {
      err <- paste0("Error: In experiment ", experiment, ", some condtions have less than 2 replicates -> ", DESeq2_condition)
      print(err)
    }, .(experiment, DESeq2_condition)]
    stop("DESeq2 needs at least two replicates to function. EXIT")
  }
  
  # Check all sample names are unique (within experiment) ----
  dups <- meta[, .(.N, count= paste0(count_table, collapse = "; ")), .(experiment, DESeq2_name)][N>1]
  if(nrow(dups))
  {
    dups[, {
      err <- paste0("Error: In experiment ", experiment, ", some sample names are associated to >1 count file: ", DESeq2_name, " -> ", count, "")
      print(err)
    }, .(experiment, DESeq2_name)]
    stop("Please update your metadata to avoid such conflict. EXIT")
  }
  
  # Print quick sum-up ----
  print(paste(length(unique(meta$experiment)), "distinct experiments were detected and will be analyzed in parallel:"))
  meta[, {
    print(paste0(experiment, " -> ", paste0(DESeq2_name, collapse = "; ")))
  }, experiment]
  
  # Create output file names for count tables and dds ----
  # Subfolder and contrast
  meta[, subFolder:= paste0(experiment, "/")]
  meta[, contrast:= paste0("_", DESeq2_condition, "_vs_", DESeq2_control)]
  # dds files
  meta[, ddsFile:= paste0(dds_output_folder, subFolder, experiment, "_DESeq2.dds")]
  # FC tables
  meta[DESeq2_condition!=DESeq2_control, fcTable:= paste0(FC_output_folder, subFolder, experiment)]
  meta[DESeq2_condition!=DESeq2_control, fcTable:= paste0(fcTable, contrast, "_DESeq2.txt")]
  # PDF MA plots
  meta[!is.na(fcTable), MAplot:= paste0(PDF_output_folder, subFolder, "MA_plots/", experiment)]
  meta[!is.na(fcTable), MAplot:= paste0(MAplot, contrast, "_DESeq2_MA_plot.pdf")]
  # Clean
  meta$subFolder <- meta$contrast <- NULL
  
  # Save processed FC metadata ----
  if(!grepl(".rds$", FC_metadata_output))
    stop("FC_metadata_output should end up with a .rds extension") else
      saveRDS(meta, FC_metadata_output)
  
  # Create directories ----
  dirs <- c(unique(dirname(na.omit(unlist(meta[, ddsFile:MAplot])))), logs)
  if(any(!dir.exists(dirs)))
  {
    sapply(dirs, dir.create, showWarnings = F, recursive = T)
    outDirs <- paste0(c(dds_output_folder, FC_output_folder, PDF_output_folder, logs), collapse= ", ")
    print(paste("Output directories were created in:", outDirs, "!"))
  }
  
  # DESeq2 commands ---- 
  meta[, DESeq2_cmd:= {
    # [required] 1/ A comma-separated list of count (ref genome)\n
    # [required] 2/ A comma-separated list of sample names \n
    # [required] 3/ A comma-separated list of conditions \n
    # [required] 4/ A comma-separated list of controls \n
    # [required] 5/ dds output folder \n
    # [required] 6/ FC tables output folder \n
    # [required] 7/ PDF output folder \n
    # [required] 8/ Experiment \n
    paste(Rpath,
          system.file("RNAseq_pipeline", "DESeq2_analysis.R", package = "vlfunctions"),
          paste0(count_table, collapse = ","),
          paste0(DESeq2_name, collapse = ","),
          paste0(DESeq2_condition, collapse = ","),
          paste0(DESeq2_control, collapse = ","),
          dds_output_folder,
          FC_output_folder,
          PDF_output_folder,
          experiment)
  }, experiment]
  
  # Return commands and submit ----
  cmd <- meta[, .(cmd= paste0(unique(na.omit(unlist(.SD))), collapse= "; ")), experiment, .SDcols= patterns("_cmd$")]
  if(nrow(cmd))
  {
    if(submit)
      cmd[, {
        vl_bsub(cmd, 
                cores= cores, 
                m = mem, 
                name = "RNAseq", 
                t = time,
                o= paste0(normalizePath(logs), "/", experiment),
                e= paste0(normalizePath(logs), "/", experiment))
      }, .(experiment, cmd)] else
        return(cmd)
  }else
    warning("No commands detected. Is the metadata empty?")
}
