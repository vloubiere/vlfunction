#' ORFeome pipeline
#' 
#' @export
vl_ORFeome_processing <- function(metadata, ...) UseMethod("vl_ORFeome_processing")

#' @describeIn vl_ORFeome_processing for excel files path
#' @export
vl_ORFeome_processing.character <- function(metadata,
                                            processed_metadata_output= gsub(".xlsx$", "_processed.rds", metadata),
                                            ...)
{
  sheet <- readxl::read_xlsx(metadata)
  start <- cumsum(sheet[[1]]=="user")>0
  sheet <- sheet[(start),]
  meta <- as.data.table(sheet[-1,])
  names(meta) <- unlist(sheet[1,])
  vl_ORFeome_processing(metadata= meta,
                        processed_metadata_output= processed_metadata_output,
                        ...)
}

#' @describeIn vl_ORFeome_processing default method
#' @export
vl_ORFeome_processing.default <- function(metadata,
                                          processed_metadata_output,
                                          fq_output_folder= "/scratch/stark/vloubiere/ORFeome/fq/",
                                          bam_output_folder= "/scratch/stark/vloubiere/ORFeome/bam/",
                                          alignment_stats_output_folder= "db/alignment_stats/ORFeome/",
                                          BC_count_output_folder= "db/counts/ORFeome/",
                                          Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                          cores= 6,
                                          mem= 32,
                                          overwrite= FALSE,
                                          counts.exists= TRUE,
                                          submit= FALSE,
                                          wdir= getwd(),
                                          logs= "db/logs/ORFeome/processing",
                                          time= '12:00:00')
{
  # Import metadata and check format ----
  meta <- data.table::copy(metadata)
  cols <- c("user", "batch", "used", "screen", "condition", "sort", "cell_line", "replicate", "sampleID", "dictionary", "MAGeCK_sort", "i7", "layout", "sequencer", "bam_path")
  if(!all(cols %in% names(meta)))
    stop(paste("Columns missing ->", paste(cols[!cols %in% names(meta)], collapse = "; ")))
  if(!all(meta[, sampleID==paste0(screen, "_", condition, "_", sort, "_", cell_line, "_", replicate)]))
    warning("sampleID should be the concatenation of screen, condition, sort, cell_line and replicate (separated by '_'). Any samples with the same sampleID will be collapsed!")
  if(!all(meta$layout %in% c("SINGLE", "PAIRED")))
    stop("layout column should only contain 'SINGLE' or 'PAIRED'")
  
  # Check used samples
  meta[, used:= as.logical(used)]
  if(any(!meta$used))
    print(paste(sum(!meta$used), "rows were set to used==FALSE and were removed"))
  meta <- meta[(used)]
  
  # Check dictionary ----
  if(any(meta$dictionary!="lib200"))
    stop("For now, only 'lib200' dictionary is supported (see 'dictionary' column of the metadata)")
  
  # Generate output paths ----
  meta[, fq1:= paste0(fq_output_folder, "/", gsub(".bam", paste0("_", i7, ifelse(layout=="PAIRED", "_1.fq.gz", ".fq.gz")), basename(bam_path))), .(bam_path, i7, layout)]
  meta[, fq1_trimmed:= gsub(".fq.gz", "_trimmed.fq.gz", fq1)]
  # re-sequencing are merged from this step on!
  meta[, bam:= paste0(bam_output_folder, "/", sampleID, "_", dictionary, ".bam")]
  meta[, align_stats:= paste0(alignment_stats_output_folder, "/", sampleID, "_", dictionary, "_stats.txt")]
  # UMI-collapsed read counts (total read counts and UMI-collapsed read counts per unique genomic coordinate)
  meta[, count:= paste0(BC_count_output_folder, "/", sampleID, "_", dictionary, "_counts.txt")]
  # Save processed metadata ----
  if(!grepl(".rds$", processed_metadata_output))
    stop("processed_metadata_output should end up with a .rds extension") else
      saveRDS(meta, processed_metadata_output)
  
  # Create output directories ----
  dirs <- c(logs,
            na.omit(unique(dirname(unlist(meta[, fq1:count])))))
  if(any(!dir.exists(dirs)))
  {
    sapply(dirs, dir.create, showWarnings = F, recursive = T)
    outDirs <- paste0(c(fq_output_folder, bam_output_folder, alignment_stats_output_folder, BC_count_output_folder), collapse= ", ")
    print(paste("Output directories were created in:", outDirs, "!"))
  }
  
  # Check whether demultiplexing and alignment should be performed ----
  meta[, preProcess:= ifelse(counts.exists, !file.exists(count), TRUE)] # Check if demultiplex/alignment should be performed
  
  # Extract fastq files ----
  meta[, extract_cmd:= {
    if(overwrite | (preProcess & !file.exists(fq1)))
    {
      fq_prefix <- normalizePath(gsub(ifelse(layout=="PAIRED", "_1.fq.gz$", ".fq.gz$"), "", fq1), mustWork = F)
      paste("samtools view -@", cores-1, bam_path,
            # "| head -n 40000", # For tests
            "| perl",
            fcase(layout=="PAIRED", system.file("ORFeome_pipeline", "demultiplex_pe.pl", package = "vlfunctions"),
                  layout=="SINGLE", system.file("ORFeome_pipeline", "demultiplex_se.pl", package = "vlfunctions")), 
            i7,
            fq_prefix)
    }
  }, .(fq1, bam_path, i7, layout, preProcess)]
  
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
  load_cmd <- paste(c(paste("cd", wdir),
                      "module load build-env/2020",
                      "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                      "module load samtools/1.9-foss-2018b",
                      "module load bowtie2/2.3.4.2-foss-2018b"), collapse = "; ")
  cols <- intersect(c("extract_cmd", "trim_cmd", "align_cmd", "counts_cmd"),
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
                name = "ORFeome", 
                t = time,
                o= paste0(normalizePath(logs), "/", sampleID),
                e= paste0(normalizePath(logs), "/", sampleID))
      }, .(sampleID, cmd)] else
        return(cmd)
  }else
    warning("All output files already existed! No command submitted ;). Consider overwrite= T if convenient.")
}

#' MAGeCK analysis ORFeome screen
#'
#' @export
vl_ORFeome_MAGeCK_auto <- function(processed_metadata, ...) UseMethod("vl_ORFeome_MAGeCK_auto")

#' @describeIn vl_ORFeome_MAGeCK_auto  for processed_metadata file path
#' @export
vl_ORFeome_MAGeCK_auto.character <- function(processed_metadata,
                                             ...)
{
  meta <- readRDS(processed_metadata)
  vl_ORFeome_MAGeCK_auto(processed_metadata= meta,
                         ...)
}

#' @describeIn vl_ORFeome_MAGeCK_auto default method
#' @export
vl_ORFeome_MAGeCK_auto.default <- function(processed_metadata,
                                           FC_metadata_output= gsub("_processed.rds", "_FC_tables.rds"),
                                           screen_group_columns= c("screen", "cell_line"),
                                           condition_group_columns= c("condition", "sort"),
                                           output_folder= "db/FC_tables/ORFeome/",
                                           input.cutoff.FUN= function(x) sum(x)>=0,
                                           sample.cutoff.FUN= function(x) sum(x)>=3,
                                           row.cutoff.FUN= function(x) sum(x)>=0,
                                           input.pseudocount= 0,
                                           sample.pseudocount= 0,
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
  # Hard copy metadata ----
  meta <- data.table::copy(processed_metadata)
  meta[, screen_name:= paste0(unlist(.BY), collapse = "_"), screen_group_columns]
  meta[, cdition_name:= paste0(unlist(.BY), collapse = "_"), condition_group_columns]
  
  # Check columns ----
  if(!"dictionary" %in% names(meta))
    stop("'dictionary' column is missing from the metadata, while required for sanity check.")
  if(!"MAGeCK_sort" %in% names(meta))
    stop("'MAGeCK_sort' column is missing from the metadata, while used to set MAGeCK 'sort' parameter ('pos' or 'neg').")
  
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
  cmb[, paired:= lengths(sample.names)==lengths(input.names)] # Samples can be paired if 
  
  # Print report ----
  message(paste(uniqueN(cmb[, screen_name]), "screen were detected and will be analyzed in parallel:"))
  cmb[, {
    message(paste0("--------------------\nScreen ", .GRP, " --> ", .N, " comparisons : "))
    .SD[, {
      .s <- paste0(unlist(sample.names), collapse = " ")
      .i <- paste0(unlist(input.names), collapse = " ")
      message(paste0(.s, " vs. ", .i))
    }, cdition_name]
  }, screen_name]
  
  # Generate output paths ----
  cmb[, raw_counts_table:= paste0(output_folder, "/", screen_name, "/", screen_name, "_", cdition_name, "_raw_counts.txt")]
  cmb[, filtered_counts_table:= gsub("_raw_counts.txt", "_filtered_counts.txt", raw_counts_table)]
  cmb[, FC_table:= gsub("_raw_counts.txt", ".gene_summary.txt", raw_counts_table)]
  cmb[, MA_plot_pdf:= gsub("_MA_plot.pdf", ".gene_summary.txt", raw_counts_table)]
  
  # Save metadata table ----
  message(paste("Metadata saved in", FC_metadata_output))
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
                      input.pseudocount= input.pseudocount,
                      sample.pseudocount= sample.pseudocount,
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

#' test
#' @export
vl_ORFeome_MAGeCK <- function(sample.counts,
                              input.counts,
                              screen_name,
                              cdition_name,
                              sample.names= gsub("^(.*?_rep\\d+).*", "\\1", basename(sample.counts)),
                              input.names= gsub("^(.*?_rep\\d+).*", "\\1", basename(input.counts)),
                              output_folder= "db/FC_tables/ORFeome/",
                              sort= "pos",
                              input.cutoff.FUN= function(x) sum(x)>=0,
                              sample.cutoff.FUN= function(x) sum(x)>=3,
                              row.cutoff.FUN= function(x) sum(x)>=0,
                              input.pseudocount= 0,
                              sample.pseudocount= 0,
                              paired= length(sample.names)==length(input.names),
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
  
  # Print report ----
  output_folder <- paste0(output_folder, "/", screen_name, "/")
  message("SCREEN ============")
  message(paste0(screen_name, ": paired=", paired, ", selection=", sort, " -> ", paste0(sample.names, collapse= ","), " vs. ", paste0(input.names, collapse= ",")))
  message(paste0("OUTPUT files saved in ", output_folder))
  message("===================\n")
  
  # Generate output files paths ----
  raw_counts_table <- paste0(output_folder, screen_name, "_", cdition_name, "_raw_counts.txt")
  filtered_counts_table <- gsub("_raw_counts.txt", "_filtered_counts.txt", raw_counts_table)
  FC_table <- gsub("_raw_counts.txt", ".gene_summary.txt", raw_counts_table)
  MA_plot_pdf <- gsub("_MA_plot.pdf", ".gene_summary.txt", raw_counts_table)
  
  # Create output directories
  dirs <- c(logs, output_folder)
  if(any(!dir.exists(dirs)))
  {
    print("Output folder were created:")
    sapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
  }
  
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
    # 8/ Pseudocount to be added to sample columns
    # 9/ Pseudocount to be added to input columns
    # 10/ Raw counts output file
    # 11/ Filtered counts output file"
    # Convert the functions to string form
    # Convert the function definitions to strings
    
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
          sample.pseudocount,
          input.pseudocount,
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
  
  # MA plot ----
  if(overwrite | !file.exists(MA_plot_pdf))
  {
    cmd <- c(cmd,
             paste(Rpath, system.file("ORFeome_pipeline", "volcano_plots_MAgECK.R", package = "vlfunctions"),
                   FC_table,
                   logFC.cutoff,
                   FDR.cutoff, 
                   MA_plot_pdf))
  }
  
  # Return commands and submit ----
  cmd <- paste0(na.omit(cmd), collapse= "; ")
  if(cmd!="")
  {
    if(submit)
    {
      vl_bsub(cmd, 
              cores= cores, 
              m = mem, 
              name = "MAGeCK", 
              t = time,
              o= paste0(normalizePath(logs), "/", screen_name, "/", cdition_name),
              e= paste0(normalizePath(logs), "/", screen_name, "/", cdition_name))
    }else
      return(cmd)
  }else
    print("All files already exist!")
}