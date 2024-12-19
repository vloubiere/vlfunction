#' ORFeome pipeline
#' 
#' To use this this pipeline, please install the vl_function package using install_github("vloubiere/vlfunction")
#' in your interactive Rstudio session AND in the local installation of R 4.3.0 ('/software/f2022/software/r/4.3.0-foss-2022b/bin/R') 
#' which will be used to submit R sub-script.
#' 
#' The pipeline is split into two main functions. The first function vl_ORFeome_processing() aligns the reads and filters confident alignments.
#' It takes as input a (correctly formatted) metadata file, saves the processed_metadata file and returns and/or submit the command lines to: \cr
#' 1/ extract reads from VBC bam file. Output .fq files are saved in fq_output_folder/ \cr
#' 2/ trim the reads. Output .fq files are saved in fq_output_folder/ \cr
#' 3/ aligns them to the library BCs (see 'library' column of the metadata table). Output .bam files are saved in bam_output_folder/ \cr
#' 4/ Compute BC counts and save the output files in BC_count_output_folder/ \cr
#' 
#' The second function, ?vl_ORFeome_MAGeCK() can be used for differential analysis. Output files are all saved in the provided output_folder/: \cr
#' 1/ Raw counts table  \cr
#' 2/ Filtered counts table  \cr
#' 3/ MAGeCK output \cr
#' 4/ MA plot .pdf \cr
#'
#' @param metadata The path to a correctly formatted .xlsx. .rds or .txt metadata file, or a data.table. See the template at '/groups/stark/vloubiere/projects/vl_pipelines/Rdata/metadata_PROseq.xlsx'.
#' @param processed_metadata_output An .rds path where to save the processed metadata file, which contains the paths of all output files and will be used to manage them.
#' By default, when importing the metadata from an excel sheet, "_processed.rds" will be appended to the excel file path. 
#' @param fq_output_folder Output folder for .fq files. Default= "/scratch/stark/vloubiere/fq/ORFeome/".
#' @param bam_output_folder Output folder for aligned bam files Default= "/scratch/stark/vloubiere/bam/ORFeome/".
#' @param alignment_stats_output_folder Output folder for alignment statistics. Default= "db/alignment_stats/ORFeome/".
#' @param BC_count_output_folder Output folder for count files. Default= "db/counts/ORFeome/".
#' @param Rpath Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of cores per job. Default= 6.
#' @param mem Memory per job (in Go). Default= 32.
#' @param overwrite Should existing files be overwritten?
#' @param counts.exists If the output umi file already exists and overwrite is set to FALSE, skips demultiplexing and alignment. Default= TRUE.
#' @param submit Should the command be submitted? Default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs/ORFeome/processing".
#' @param time The time required for the SLURM scheduler. Default= '12:00:00'.
#'
#' @return Return a data.table containing, for each sampleID, the concatenated commands that are required to process the data.
#' These commands can then be submitted either directly via the function, or using vl_bsub()...
#'
#' @examples
#' library(vlfunctions)
#' 
#' # Example metadata ----
#' metadata <- system.file("ORFeome_pipeline", "metadata_ORFeome.xlsx", package = "vlfunctions")
#' 
#' # Generate command lines ----
#' cmd <- vl_ORFeome_processing(metadata= metadata,
#'                              overwrite= TRUE,
#'                              submit= FALSE)
#' # To actually submit the jobs, switch to submit= TRUE
#' vl_ORFeome_processing(metadata= metadata,
#'                       processed_metadata_output = "Rdata/metadata_ORFeome_processed.rds",
#'                       overwrite= FALSE,  # Existing files will not be processed again
#'                       submit= TRUE)
#' 
#' # Differential analysis
#' See ?vl_ORFeome_MAGeCK()
#' 
#' @export
vl_ORFeome_processing <- function(metadata,
                                  processed_metadata_output,
                                  fq_output_folder= "/scratch/stark/vloubiere/fq/ORFeome/",
                                  bam_output_folder= "/scratch/stark/vloubiere/bam/ORFeome/",
                                  alignment_stats_output_folder= "db/alignment_stats/ORFeome/",
                                  BC_count_output_folder= "db/counts/ORFeome/",
                                  Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                  cores= 8,
                                  mem= 64,
                                  overwrite= FALSE,
                                  counts.exists= TRUE,
                                  submit= FALSE,
                                  wdir= getwd(),
                                  logs= "db/logs/ORFeome/processing",
                                  time= '12:00:00')
{
  # Import metadata and check layout ----
  meta <- if(is.data.table(metadata))
    data.table::copy(metadata) else
      vl_import_metadata_sheet(metadata)
  
  # Check required columns ----
  # Design
  required_cols <- c("used", "screen", "condition", "sort", "cell_line", "replicate",
                     "sampleID", "dictionary", "MAGeCK_sort", "layout", "fq1", "fq2")
  missing_cols <- setdiff(required_cols, names(meta))
  if(length(missing_cols))
    stop(paste("Metadata required columns mising:", paste0(missing_cols, collapse = ", ")))
  # Check fq files
  meta[, fq1:= as.character(fq1)]
  meta[, fq2:= as.character(fq2)]
  if(nrow(meta[is.na(fq1)]) | nrow(meta[layout=="PAIRED" & is.na(fq2)]))
  {
    missing_cols <- setdiff(c("bam_path", "i7"), names(meta))
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
  # Dictionary alignment indexes
  if(any(meta$dictionary!="lib200"))
    stop("For now, only 'lib200' dictionary is supported (see 'dictionary' column of the metadata)")
  # layout
  if(!all(meta$layout %in% c("SINGLE", "PAIRED")))
    stop("layout column should only contain 'SINGLE' or 'PAIRED'")
  # Check used samples
  meta[, used:= as.logical(used)]
  if(any(!meta$used))
    print(paste(sum(!meta$used), "rows were set to used==FALSE and were removed"))
  meta <- meta[(used)]
  # Warning for duplicated sample IDs
  potential_dup <- meta[, .N, sampleID][N>1]
  if(nrow(potential_dup))
    warning(paste("The following sampleIDs were present more than once and will be merged (resequenced)?:",
                  paste0(potential_dup$sampleID, collapse = ", ")))
  
  # Generate output paths ----
  # fq files
  meta[is.na(fq1) & !is.na(bam_path), fq1:= {
    paste0(fq_output_folder, "/",
           gsub(".bam", paste0("_", i7, "_1.fq.gz"), basename(bam_path)))
  }, .(bam_path, i7, layout)]
  # Trimmed
  meta[, fq1_trimmed:= paste0(fq_output_folder, "/", gsub(".fq.gz", "_trimmed.fq.gz", basename(fq1)))]
  # re-sequencing are merged from this step on!
  meta[, bam:= paste0(bam_output_folder, "/", sampleID, "_", dictionary, ".bam")]
  meta[, align_stats:= paste0(alignment_stats_output_folder, "/", sampleID, "_", dictionary, "_stats.txt")]
  # UMI-collapsed read counts (total read counts and UMI-collapsed read counts per unique genomic coordinate)
  meta[, count:= paste0(BC_count_output_folder, "/", sampleID, "_", dictionary, "_counts.txt")]
  
  # Check whether demultiplexing and alignment should be performed ----
  meta[, preProcess:= ifelse(counts.exists, !file.exists(count), TRUE)] # Check if demultiplex/alignment should be performed
  
  # Extract fastq files ----
  meta[, extract_cmd:= {
    if(overwrite | (preProcess & !file.exists(fq1)))
    {
      fq_prefix <- normalizePath(gsub("_1.fq.gz$", "", fq1), mustWork = F)
      paste("samtools view -@", cores-1, bam_path,
            # "| head -n 40000", # For tests
            "| perl",
            fcase(layout=="PAIRED", system.file("ORFeome_pipeline", "demultiplex_pe.pl", package = "vlfunctions"),
                  layout=="SINGLE", system.file("ORFeome_pipeline", "demultiplex_se.pl", package = "vlfunctions")), 
            i7,
            fq_prefix)
    }
  }, .(layout, fq1, preProcess)]
  
  # Trim fq files ----
  meta[, trim_cmd:= {
    if(overwrite | (preProcess & !file.exists(fq1_trimmed)))
    {
      paste0("trim_galore --gzip -o ", dirname(fq1_trimmed), "/ ", fq1)
    }
  }, .(fq1, fq1_trimmed, preProcess)]
  
  # Align to BC sequences ----
  meta[, align_cmd:= {
    if(overwrite | (preProcess & !file.exists(bam)))
    {
      # Bowtie index BARCODES
      idx <- switch(dictionary,
                    "lib200"= "/groups/stark/vloubiere/projects/viralORF_tomas/db/indexes_BCs/lib200_merged/ORF")
      
      # Command
      paste("bowtie2 -p 6 -U",
            paste0(fq1_trimmed, collapse= ","),
            "-x", idx, 
            "| samtools sort -@", cores-1, "-o", bam)
    }
  }, .(bam, dictionary, preProcess)]
  
  # UMI counts ----
  meta[, counts_cmd:= {
    if(overwrite | !file.exists(count))
    {
      # Retrieve BC .rds file ----
      BC <- switch(dictionary,
                   "lib200"= "/groups/stark/vloubiere/projects/viralORF_tomas/db/dictionary/lib200_merged_dictionary.rds")
      
      # Compute BC counts ----
      paste(Rpath,
            system.file("ORFeome_pipeline", "BC_counts.R", package = "vlfunctions"),
            bam,
            BC,
            align_stats,
            count)
    }
  }, .(bam, dictionary, align_stats, count)]
  
  # Return commands ----
  load_cmd <- paste(c("module load build-env/2020",
                      "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                      "module load samtools/1.9-foss-2018b",
                      "module load bowtie2/2.3.4.2-foss-2018b"), collapse = "; ")
  cols <- intersect(c("extract_cmd", "trim_cmd", "align_cmd", "counts_cmd"),
                    names(meta))
  cmd <- if(length(cols))
    meta[, .(cmd= paste0(c(load_cmd, unique(na.omit(unlist(.SD)))),
                         collapse = "; ")), sampleID, .SDcols= cols][cmd!=load_cmd] else
                           data.table()
  
  # If commands were generated ----
  if(nrow(cmd))
  {
    # If commands are to be submitted ----
    if(submit)
    {
      # Save processed metadata ----
      saveRDS(meta[, -c(cols, "preProcess"), with= FALSE],
              processed_metadata_output)
      
      # Create output directories ----
      dirs <- c(logs,
                na.omit(unique(dirname(unlist(meta[, fq1:count])))))
      if(any(!dir.exists(dirs)))
      {
        sapply(dirs, dir.create, showWarnings = F, recursive = T)
        outDirs <- paste0(c(fq_output_folder,
                            bam_output_folder,
                            alignment_stats_output_folder,
                            BC_count_output_folder),
                          collapse= ", ")
        print(paste("Output directories were created in:", outDirs, "!"))
      }
      
      # Submit commands ----
      cmd[, {
        vl_bsub(cmd, 
                cores= cores, 
                m = mem, 
                name = "ORFeome", 
                t = time,
                o= paste0(normalizePath(logs), "/", sampleID),
                e= paste0(normalizePath(logs), "/", sampleID),
                wdir= wdir)
      }, .(sampleID, cmd)] 
    }else
    {
      return(cmd)
    }
  }else
    warning("All output files already existed! No command submitted ;). Consider overwrite= T if convenient.")
}

#' MAGeCK analysis ORFeome screen
#'
#' This function uses the provided metadata file to run MAGeCK comparisons, which can be automatized using ?vl_ORFeome_MAGeCK_auto().
#'
#' @param sample.counts Vector of sample count files.
#' @param input.counts Vector of input count files.
#' @param screen_name Screen name (for example 'ISRE_A549'), which will be used to create a subfolder in output_folder.
#' @param cdition_name Condition name (for example 'IFN_dim'), which will be used to create sample output name.
#' @param sample.names Sample names, for example c("sort_rep1", "sort_rep2").
#' @param input.names Sample names, for example c("input_rep1", "input_rep2").
#' @param output_folder Output folder where all output files should be saved. Default= "db/FC_tables/ORFeome/"
#' @param sort Sort parameter of the screen, either 'pos' or 'neg'. Default= 'pos'.
#' @param sample.cutoff.FUN Function to be applied to filter sample columns. Default= function(x) sum(x)>=3
#' @param input.cutoff.FUN Function to be applied to filter input columns. Default= function(x) sum(x)>=0
#' @param row.cutoff.FUN Function to be applied to filter all columns. Default= function(x) sum(x)>=3
#' @param pseudocount The pseudocount to be added to 0 values (only!) in sample and input count columns (note that this value will be further normalized for sequencing depth). Default= 0 
#' @param paired Should samples be analyzed as paired? Default= FALSE.
#' @param logFC.cutoff logFC cutoff to call hits for the MA plot. Default= 1
#' @param FDR.cutoff FDR cutoff to call hits for the MA plot. Default= 0.05
#' @param Rpath Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of cores per job. Default= 2.
#' @param mem Memory per job (in Go). Default= 8.
#' @param overwrite Should existing files be overwritten?
#' @param submit Should the command be submitted? Default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs/ORFeome/MAGeCK".
#' @param time The time required for the SLURM scheduler. Default= '1:00:00'.
#'
#' @return Command lines to run MAGeCK for a single screen.
#' @export
vl_ORFeome_MAGeCK <- function(sample.counts,
                              input.counts,
                              screen_name,
                              cdition_name,
                              sample.names= gsub("^(.*?_rep\\d+).*", "\\1", basename(sample.counts)),
                              input.names= gsub("^(.*?_rep\\d+).*", "\\1", basename(input.counts)),
                              output_folder= "db/FC_tables/ORFeome/",
                              sort= "pos",
                              sample.cutoff.FUN= function(x) sum(x)>=3,
                              input.cutoff.FUN= function(x) sum(x)>=0,
                              row.cutoff.FUN= function(x) sum(x)>=3,
                              pseudocount= 0,
                              paired= FALSE,
                              logFC.cutoff= 1,
                              FDR.cutoff= 0.05,
                              Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                              cores= 2,
                              mem= 8,
                              overwrite= FALSE,
                              submit= FALSE,
                              wdir= getwd(),
                              logs= "db/logs/ORFeome/MAGeCK",
                              time= '01:00:00')
{
  # Checks ----
  if(any(duplicated(sample.names)))
    stop("Some sample.names are duplicated. Make sure that sample.counts files are unique and provide unique sample.names.")
  if(any(duplicated(input.names)))
    stop("Some input.names are duplicated. Make sure that input.counts files are unique and provide unique sample.names.")
  if(!sort %in% c("pos", "neg"))
    stop("sort should either be set to 'pos' or 'neg'.")
  
  # Print report ----
  output_folder <- paste0(output_folder, "/", screen_name, "/")
  message("SCREEN ============")
  message(paste0(screen_name, ": paired=", paired, ", selection=", sort, " -> ", paste0(sample.names, collapse= ","), " vs. ", paste0(input.names, collapse= ",")))
  message(paste0("OUTPUT files saved in ", output_folder))
  message("===================\n")
  
  # Generate output files paths ----
  raw_counts_table <- paste0(output_folder, screen_name, "_", cdition_name, ".", sort, ".raw_counts.txt")
  filtered_counts_table <- gsub(".raw_counts.txt$", ".filtered_counts.txt", raw_counts_table)
  FC_table <- gsub(".raw_counts.txt$", ".gene_summary.txt", raw_counts_table)
  MA_plot_pdf <- gsub(".raw_counts.txt$", ".MA_plot.pdf", raw_counts_table)
  master_table_prefix <- gsub("gene_summary.txt$", "gene_summary_master", FC_table)
  
  # Generate raw and filtered counts table ----
  cmd <- if(overwrite | !file.exists(raw_counts_table) | !file.exists(filtered_counts_table))
  {
    # 1/ A comma-separated list of sample count files
    # 2/ A comma-separated list of input count files
    # 3/ A comma-separated list of sample names (mathcing the list of count files)
    # 4/ A comma-separated list of input names (mathcing the list of count files)
    # 5/ Function to be applied to input columns for filtering
    # 6/ Function to be applied to sample columns for filtering
    # 7/ Function to be applied to all columns for filtering
    # 8/ Pseudocount to be added to sample and input columns (will be further normalized for sequencing depth)
    # 9/ Raw counts output file
    # 10/ Filtered counts output file
    
    # Deparse functions
    sample_fun_string <- paste(deparse(sample.cutoff.FUN), collapse = " ")
    input_fun_string <- paste(deparse(input.cutoff.FUN), collapse = " ")
    row_fun_string <- paste(deparse(row.cutoff.FUN), collapse = " ")
    
    # Construct the system call
    paste(Rpath,
          system.file("ORFeome_pipeline", "compute_count_tables.R", package = "vlfunctions"),
          paste0(sample.counts, collapse = ","),
          paste0(input.counts, collapse = ","),
          paste0(sample.names, collapse = ","),
          paste0(input.names, collapse = ","),
          shQuote(sample_fun_string),
          shQuote(input_fun_string),
          shQuote(row_fun_string),
          pseudocount,
          raw_counts_table,
          filtered_counts_table)
  }else
    as.character(NA)
  
  # MAGeCK analysis ----
  if(overwrite | !file.exists(FC_table))
  {
    cmd <- c(cmd, 
             paste("module load build-env/2020; module load mageck/0.5.9-foss-2018b; mageck test",
                   "-k", filtered_counts_table,
                   "-t", paste0(unlist(sample.names), collapse = ","),
                   "-c", paste0(unlist(input.names), collapse = ","),
                   ifelse(paired, "--paired", ""),
                   "-n", gsub(".gene_summary.txt", "", FC_table),
                   "--sort-criteria", sort,
                   "--remove-zero none"))
  }
  
  # MA plot and merged table ----
  if(overwrite | !file.exists(MA_plot_pdf))
  {
    cmd <- c(cmd,
             paste(Rpath,
                   system.file("ORFeome_pipeline", "volcano_plots_MAgECK.R", package = "vlfunctions"),
                   FC_table,
                   logFC.cutoff,
                   FDR.cutoff, 
                   MA_plot_pdf))
  }
  
  # Merge master table plot (always overwrite cause very fast) ----
  cmd <- c(cmd,
           paste(Rpath,
                 system.file("ORFeome_pipeline", "merge_gene_summary_to_master_table.R", package = "vlfunctions"),
                 FC_table,
                 master_table_prefix,
                 sort))
  
  # Return commands and submit ----
  cmd <- paste0(na.omit(cmd), collapse= "; ")
  
  # If command had been generated ----
  if(cmd!="")
  {
    # If command is to be submitted ----
    if(submit)
    {
      # Create output directories
      dirs <- c(logs, output_folder)
      if(any(!dir.exists(dirs)))
      {
        print("Output folder were created:")
        sapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
      }
      
      # Submit commands ----
      vl_bsub(cmd, 
              cores= cores, 
              m = mem, 
              name = "MAGeCK", 
              t = time,
              o= paste0(normalizePath(logs), "/", screen_name, "/", cdition_name),
              e= paste0(normalizePath(logs), "/", screen_name, "/", cdition_name),
              wdir= wdir)
    }else
    {
      return(cmd)
    }
  }else
    print("Files already existed! No command generated.")
}

#' Automatic wrapper MAGeCK analysis ORFeome screen
#'
#' This function is a wrapper around ?vl_ORFeome_MAGeCK() that used the provided metadata to run all comparisons.
#' @param processed_metadata Path to the metadata file generated by vl_ORFeome_processing (in .rds or .txt format), or the corresponding data.table.
#' @param FC_metadata_output An .rds path where to save the FC metadata file, which contains the path to MAGeCK output files. By default, when the processed_metadata is a path to a processed_metadata files, "_FC_tables.rds" will be appended processed_metadata file path. 
#' @param screen_group_columns Columns to be used to identify individual screens. Default= c("screen", "cell_line"),
#' @param condition_group_columns Columns to be used to identify individual conditions. Default= c("screen", "cell_line"),
#' @param paired Should samples be analyzee as paired? Possible value are "auto" (which will consider samples paired if the number of sample count files matches the number of input count files), TRUE or FALSE. Default= 'auto'.
#' @param sort_column The name of the column containing sort options ('pos' or 'neg'). Default= 'MAGeCK_sort'.
#' @param output_folder Output folder where all output files should be saved. Default= "db/FC_tables/ORFeome/"
#' @param sample.cutoff.FUN Function to be applied to filter sample columns. Default= function(x) sum(x)>=3
#' @param input.cutoff.FUN Function to be applied to filter input columns. Default= function(x) sum(x)>=0
#' @param row.cutoff.FUN Function to be applied to filter all columns. Default= function(x) sum(x)>=3
#' @param pseudocount The pseudocount to be added to 0 values (only!) in sample and input count columns (note that this value will be further normalized for sequencing depth). Default= 0 
#' @param logFC.cutoff logFC cutoff to call hits for the MA plot. Default= 1
#' @param FDR.cutoff FDR cutoff to call hits for the MA plot. Default= 0.05
#' @param Rpath Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of cores per job. Default= 2.
#' @param mem Memory per job (in Go). Default= 8.
#' @param overwrite Should existing files be overwritten?
#' @param submit Should the command be submitted? Default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs/ORFeome/MAGeCK".
#' @param time The time required for the SLURM scheduler. Default= '1:00:00'.
#'
#' @return Command lines to run MAGeCK for a single screen.
#' @examples
#' # Processed metadata ----
#' processed <- readRDS("Rdata/metadata_ORFeome_processed.rds")
#' 
#' # Call hits for all conditions using MAgECK ----
#' vl_ORFeome_MAGeCK_auto(processed,
#'                        #' FC_metadata_output = "Rdata/metadata_ORFeome_FC_tables.rds",
#'                        overwrite = FALSE, # Do not re-process existing files
#'                        submit = TRUE)
#'                        
#' # FC tables metadata ----
#' FC_tables <- readRDS("Rdata/metadata_ORFeome_FC_tables.rds")
#' 
#' @export
vl_ORFeome_MAGeCK_auto <- function(processed_metadata,
                                   FC_metadata_output= gsub("_processed.rds", "_FC_tables.rds"),
                                   screen_group_columns= c("screen", "cell_line"),
                                   condition_group_columns= c("condition", "sort"),
                                   paired= "auto",
                                   sort_column= "MAGeCK_sort",
                                   output_folder= "db/FC_tables/ORFeome/",
                                   input.cutoff.FUN= function(x) sum(x)>=0,
                                   sample.cutoff.FUN= function(x) sum(x)>=3,
                                   row.cutoff.FUN= function(x) sum(x)>=0,
                                   pseudocount= 0,
                                   logFC.cutoff= 1,
                                   FDR.cutoff= 0.05,
                                   Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                   cores= 2,
                                   mem= 8,
                                   overwrite= FALSE,
                                   submit= FALSE,
                                   wdir= getwd(),
                                   logs= "db/logs/ORFeome/MAGeCK",
                                   time= '01:00:00')
{
  # Import metadata and check layout ----
  meta <- if(is.data.table(processed_metadata))
    data.table::copy(processed_metadata) else
      vl_import_metadata_sheet(processed_metadata)
  
  # screen_name, cdition_name ----
  meta[, screen_name:= paste0(unlist(.BY), collapse = "_"), screen_group_columns]
  meta[, cdition_name:= paste0(unlist(.BY), collapse = "_"), condition_group_columns]
  
  # Check columns ----
  if(!"dictionary" %in% names(meta))
    stop("'dictionary' column is missing from the metadata, while required for sanity check.")
  if(sort_column!="MAGeCK_sort")
  {
    message(paste0("MAGeCK_sort column was replaced by ", sort_column))
    meta$MAGeCK_sort <- meta[[sort_column]]
  }
  if(!all(meta[!grepl("input", cdition_name), MAGeCK_sort] %in% c("pos", "neg")))
    stop("For sorted (non-input) samples, 'MAGeCK_sort' column should either be set to 'pos' or 'neg'.")
  
  # Check metadata output
  if(submit && !grepl(".rds$", FC_metadata_output))
    stop("FC_metadata_output should end up with a .rds extension")
  
  # Make sure that replicates will be in the same order ----
  setorderv(meta, "replicate")
  
  # Collapse re-sequencing ----
  meta <- unique(meta[, .(screen_name, cdition_name, sampleID, count, MAGeCK_sort)])
  
  # Check comparisons ----
  .sample <- meta[!grepl("input", cdition_name), .(ID= .(sampleID), count= .(count)), .(screen_name, cdition_name, MAGeCK_sort)]
  setnames(.sample, c("ID", "count"), c("sample.names", "sample.counts"))
  .input <- meta[grepl("input", cdition_name), .(ID= .(sampleID), count= .(count)), .(screen_name, cdition_name)]
  setnames(.input, c("ID", "count"), c("input.names", "input.counts"))
  cmb <- merge(x = .sample,
               y = .input[, !"cdition_name"],
               by= "screen_name")
  
  # Check if analyses can be paired ----
  if(paired=="auto")
    cmb[, paired:= lengths(sample.names)==lengths(input.names)] else if(is.logical(paired))
      cmb$paired <- paired else
        stop("paired should be set to 'auto' or a logical value")
  
  # Print report ----
  message(paste(uniqueN(cmb[, screen_name]), "screen were detected and will be analyzed in parallel:"))
  cmb[, {
    message(paste0("--------------------\nScreen ", .GRP, " --> ", .N, " comparisons : "))
    .SD[, {
      .s <- paste0(unlist(sample.names), collapse = " ")
      .i <- paste0(unlist(input.names), collapse = " ")
      message(paste0("Sort: ", MAGeCK_sort, " -> ", .s, " vs. ", .i))
    }, .(MAGeCK_sort, cdition_name)]
  }, screen_name]
  
  # Generate output paths ----
  cmb[, raw_counts_table:= paste0(output_folder, "/", screen_name, "/", screen_name, "_", cdition_name, ".", MAGeCK_sort, ".raw_counts.txt"), MAGeCK_sort]
  cmb[, filtered_counts_table:= gsub(".raw_counts.txt$", ".filtered_counts.txt", raw_counts_table)]
  cmb[, FC_table:= gsub(".raw_counts.txt$", ".gene_summary.txt", raw_counts_table)]
  cmb[, MA_plot_pdf:= gsub(".raw_counts.txt$", ".MA_plot.pdf", raw_counts_table)]
  cmb[, master_table_prefix:= gsub("gene_summary.txt$", "gene_summary_master", FC_table)]
  
  # Save metadata table ----
  saveRDS(cmb,
          FC_metadata_output)
  
  # Compute comparisons ----
  cmb[, {
    .sc <- unlist(sample.counts)
    .ic <- unlist(input.counts)
    .sn <- unlist(sample.names)
    .in <- unlist(input.names)
    vl_ORFeome_MAGeCK(sample.counts= .sc,
                      input.counts=  .ic,
                      sample.names=  .sn,
                      input.names=   .in,
                      screen_name= screen_name,
                      cdition_name= cdition_name,
                      output_folder= output_folder,
                      sort= MAGeCK_sort,
                      input.cutoff.FUN= input.cutoff.FUN,
                      sample.cutoff.FUN= sample.cutoff.FUN,
                      row.cutoff.FUN= row.cutoff.FUN,
                      pseudocount= pseudocount,
                      paired= paired,
                      logFC.cutoff= logFC.cutoff,
                      FDR.cutoff= FDR.cutoff,
                      Rpath= Rpath,
                      cores= cores,
                      mem= mem,
                      overwrite= overwrite,
                      submit= submit,
                      wdir= wdir,
                      logs= logs,
                      time= time)
  }, .(screen_name, cdition_name, MAGeCK_sort, paired)]
}

