#' PRO-seq pipeline
#' 
#' To use this this pipeline, please install the vl_function package using install_github("vloubiere/vlfunction")
#' in your interactive Rstudio session AND in the local installation of R 4.3.0 ('/software/f2022/software/r/4.3.0-foss-2022b/bin/R') 
#' which will be used to submit R sub-script.
#' 
#' The pipeline is split into two main function The first function vl_PROseq_processing() aligns the reads and filters confident alignments.
#' It takes as input a (correctly formatted) metadata file, saves the processed metadata file and returns and/or submit the command lines to: \cr
#' 1/ extract reads from VBC bam file \cr
#' 2/ trim the reads \cr
#' 3/ aligns them to mouse/human genome (see 'genome' column of the metadata table) \cr
#' 4/ store unaligned read into a fastq file \cr
#' 5/ align these remaining reads to spike-in genome ('dm3') \cr
#' 6/ compute counts and statistics both for the reference genome and for spike-in \cr
#' 7/ count the reads overlapping promoter, gene-body and transcripts \cr
#' 8/ output .bw tracks for positive and negative strand reads.
#' 
#' The second function, vl_PROseq_DESeq2() can be used for differential analysis. It returns: \cr
#' 1/ a .dds DESeq2 object for each experiment \cr
#' 2/ the fold-change tables corresponding to the comparison of 'DESeq2_condition' versus 'DESeq2_control' conditions (see colnames of the processed_metadata file) \cr
#' 3/ A barplot showing the statistics of each dds file \cr
#' 4/ MA plots corresponding to each comparison \cr
#'
#' @param metadata The path to a correctly formatted .xlsx metadata file or a data.table. See the template at '/groups/stark/vloubiere/projects/vl_pipelines/Rdata/metadata_PROseq.xlsx'.
#' @param processed_metadata_output An .rds path where to save the processed metadata file, which contains the paths of all output files and will be used to manage them.
#' By default, when importing the metadata from an excel sheet, "_processed.rds" will be appended to the excel file path. 
#' @param count_output_folder Output folder for count files. Default= "db/counts/PROseq/".
#' @param alignment_stats_output_folder Output folder for alignment statistics. Default= "db/alignment_stats/PROseq/".
#' @param count_table_output_folder Output folder for count tables. Default= "db/count_tables/PROseq/".
#' @param bw_output_folder Output folder for .bw tracks. Default= "db/bw/PROseq/".
#' @param tmp_folder Output folder for temporary files (.fq, .bam). Default= "/scratch/stark/vloubiere".
#' @param Rpath Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param cores Number of cores per job. Default= 8.
#' @param mem Memory per job (in Go). Default= 64.
#' @param overwrite Should existing files be overwritten?
#' @param submit Should the command be submitted? Default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs/CUTNRUN/processing".
#' @param time The time required for the SLURM scheduler. Default= '1-00:00:00'.
#'
#' @return Return a data.table containing, for each sampleID, the concatenated commands that are required to process the data.
#' These commands can then be submitted either directly via the function, or using vl_bsub()...
#'
#' @examples
#' # Process example dataset ----
#' meta <- vl_metadata_PROseq
#' vl_PROseq_processing(metadata = meta,
#'                      processed_metadata_output = "Rdata/metadata_PROseq_processed.rds",
#'                      cores= 8,
#'                      mem= 64,
#'                      overwrite= FALSE,
#'                      submit= TRUE)
#' 
#' # Once the jobs have finished to work, you can check the size of output files:
#' processed <- readRDS("Rdata/metadata_PROseq_processed.rds")
#' file.size(na.omit(unlist(processed[, fq1:bwNS])))
#' 
#' # Differential analysis (of note, no "overwrite" option is available in this case)
#' vl_PROseq_DESeq2(processed, # Here, you could also provide the path to the processed metadata table ("Rdata/metadata_PROseq_processed.rds").
#'                  FC_metadata_output = "Rdata/metadata_PROseq_FC_tables.rds",
#'                  submit = TRUE)
#'                  
#' # Once the commands ran, to retrieve spikein normalized FC tables:
#' FC <- readRDS("Rdata/metadata_PROseq_FC_tables.rds")
#' files <- FC[, list.files(fcDir, paste0(DESeq2_condition, ".*", "spike")), .(fcDir, DESeq2_condition)]$V1
#'
#' @export
vl_PROseq_processing <- function(metadata, ...) UseMethod("vl_PROseq_processing")

#' @describeIn vl_PROseq_processing for excel files path
#' @export
vl_PROseq_processing.character <- function(metadata,
                                           processed_metadata_output= gsub(".xlsx$", "_processed.rds", metadata),
                                           ...)
{
  sheet <- readxl::read_xlsx(metadata)
  start <- cumsum(sheet[[1]]=="user")>0
  sheet <- sheet[(start),]
  meta <- as.data.table(sheet[-1,])
  names(meta) <- unlist(sheet[1,])
  vl_PROseq_processing(metadata= meta,
                       processed_metadata_output= processed_metadata_output,
                       ...)
}

#' @describeIn vl_PROseq_processing default method
#' @export
vl_PROseq_processing.default <- function(metadata,
                                         processed_metadata_output,
                                         count_output_folder= "db/counts/PROseq/",
                                         alignment_stats_output_folder= "db/alignment_stats/PROseq/",
                                         count_table_output_folder= "db/count_tables/PROseq/",
                                         bw_output_folder= "db/bw/PROseq/",
                                         tmp_folder= "/scratch/stark/vloubiere",
                                         Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                         cores= 8,
                                         mem= 64,
                                         overwrite= FALSE,
                                         submit= FALSE,
                                         wdir= getwd(),
                                         logs= "db/logs/PROseq/processing",
                                         time= '1-00:00:00')
{
  # Import metadata and check format ----
  meta <- data.table::copy(metadata)
  cols <- c("user", "batch", "experiment", "cell_line", "treatment", "replicate", "sampleID",
            "DESeq2_name", "DESeq2_condition", "DESeq2_control", "barcode", "eBC", "genome",
            "spikein_genome", "layout", "sequencer", "bam_path")
  if(!all(cols %in% names(meta)))
    stop(paste("Columns missing ->", paste(cols[!cols %in% names(meta)], collapse = "; ")))
  if(any(!na.omit(meta$DESeq2_control) %in% meta$DESeq2_condition))
    stop("Some DESeq2_control values have no correspondence in the DESeq2_condition column!")
  if(!all(meta[, sampleID==paste0(cell_line, "_", treatment, "_", experiment, "_", replicate)]))
    warning("sampleID should be the concatenation of cell_line, treatment, experiment and replicate (separated by '_'). Any samples with the same sampleID will be collapsed!")
  if(!all(meta$layout %in% c("SINGLE", "PAIRED")))
    stop("layout column should only contain 'SINGLE' or 'PAIRED'")
  
  # Check bowtie2 idx ----
  if(any(!meta$genome %in% c("mm10", "hg38")))
    stop("Only mm10 and hg38 are supported! For other genomes, please provide path to the corresponding bowtie 2 index.")
  if(any(!meta$spikein_genome %in% c("dm3")))
    stop("Only dm3 is supported for spike-in! For other genomes, please provide path to the corresponding bowtie 2 index.")
  
  # Generate output paths ----
  meta[, fq1:= paste0(tmp_folder, "/PROseq/fq/", gsub(".bam", paste0("_", barcode, "_", eBC, "_1.fq.gz"), basename(bam_path))), .(bam_path, barcode, eBC)]
  meta[, fq1_trimmed:= gsub(".fq.gz", "_trimmed.fq.gz", fq1)]
  # re-sequencing are merged from this step on!
  meta[, bam:= paste0(tmp_folder, "/PROseq/bam/", sampleID, "_", genome, "_aligned.bam")]
  # Unaligned reads alignment to spike-in genome
  meta[, fq_unaligned:= paste0(tmp_folder, "/PROseq/fq/", sampleID, "_unaligned.fq")] 
  meta[, bamSpike:= paste0(tmp_folder, "/PROseq/bam/", sampleID, "_", spikein_genome, "_spikein_aligned.bam")]
  # UMI-collapsed read counts (total read counts and UMI-collapsed read counts per unique genomic coordinate)
  meta[, count:= paste0(count_output_folder, experiment, "/", sampleID, "_", genome, "_counts.txt")]
  meta[, countSpike:= paste0(count_output_folder, experiment, "/", sampleID, "_", spikein_genome, "_spikein_counts.txt")]
  # reads statistics
  meta[, read_stats:= paste0(alignment_stats_output_folder, experiment, "/", sampleID, "_", genome, "_statistics.txt")]
  meta[, spike_stats:= paste0(alignment_stats_output_folder, experiment, "/", sampleID, "_", spikein_genome, "_spikeIn_statistics.txt")]
  # Aggregate read counts per gene feature (TSS, body, transcript)
  meta[, countTables_promoter:= paste0(count_table_output_folder, experiment, "/promoter/", sampleID, "_", genome, "_counts.txt")]
  meta[, countTables_geneBody:= paste0(count_table_output_folder, experiment, "/gene_body/", sampleID, "_", genome, "_counts.txt")]
  meta[, countTables_transcript:= paste0(count_table_output_folder, experiment, "/transcript/", sampleID, "_", genome, "_counts.txt")]
  # bw tracks
  meta[, bwPS:= paste0(bw_output_folder, experiment, "/", DESeq2_name, ".ps.bw")]
  meta[, bwNS:= paste0(bw_output_folder, experiment, "/", DESeq2_name, ".ns.bw")]
  
  # Save processed metadata ----
  if(!grepl(".rds$", processed_metadata_output))
    stop("processed_metadata_output should end up with a .rds extension") else
      saveRDS(meta, processed_metadata_output)
  
  # Create output directories ----
  dirs <- c(logs,
            na.omit(unique(dirname(unlist(meta[, fq1:bwNS])))))
  if(any(!dir.exists(dirs)))
  {
    sapply(dirs, dir.create, showWarnings = F, recursive = T)
    outDirs <- paste0(c(count_output_folder, count_table_output_folder, bw_output_folder, tmp_folder), collapse= ", ")
    print(paste("Output directories were created in:", outDirs, "!"))
  }
  
  # Extract fastq files ----
  meta[, extract_cmd:= {
    if(overwrite | !file.exists(fq1))
    {
      fq_prefix <- normalizePath(gsub("_1.fq.gz$", "", fq1), mustWork = F)
      paste("samtools view -@", cores-1, bam_path,
            # "| head -n 40000", # For tests
            "| perl",
            fcase(layout=="PAIRED", system.file("PROseq_pipeline", "PROseq_demultiplexing_pe.pl", package = "vlfunctions"), # eBC length= 4; UMI length= 10
                  layout=="SINGLE", system.file("PROseq_pipeline", "PROseq_demultiplexing_se.pl", package = "vlfunctions")), 
            barcode,
            eBC,
            fq_prefix, "; gzip", paste0(fq_prefix, "_1.fq"))
    }
  }, .(fq1, bam_path, barcode, eBC, layout)]
  
  # Trim fq files ----
  meta[, trim_cmd:= {
    if(overwrite | !file.exists(fq1_trimmed))
    {
      paste("cutadapt -a TGGAATTCTCGGGTGCCAAGG",
            "-m", 10, # Remove trimmed reads shorter than 10bp
            "-o", fq1_trimmed, fq1)
    }
  }, .(fq1, fq1_trimmed)]
  
  # Align to reference genome ----
  meta[, align_cmd:= {
    if(overwrite | !file.exists(bam))
    {
      # Bowtie index reference genome
      refIdx <- switch(genome,
                       "mm10"= "/groups/stark/indices/bowtie/mm10/mm10",
                       "hg38"= "/groups/stark/indices/bowtie/hg38/hg38")
      # Command
      paste("bowtie -q -v 2 -m 1 --best --strata --sam -p", cores,
            refIdx,
            paste0(fq1_trimmed, collapse = ","),
            "| samtools view -bS -o ", bam)
    }
  }, .(bam, genome)]
  
  # Extract unmapped reads ----
  meta[, unaligned_cmd:= {
    if(overwrite | !file.exists(fq_unaligned))
      paste("samtools view -f 4 -b", bam, "| samtools fastq - >", fq_unaligned)
  }, .(bam, fq_unaligned)]
  
  # Align unmapped reads to spike-in genome ----
  meta[, alignSpike_cmd:= {
    if(overwrite | !file.exists(bamSpike))
    {
      # Bowtie index spike in genome
      spikeIdx <- switch(spikein_genome,
                         "dm3"= "/groups/stark/indices/bowtie/dm3/dm3")
      # Command
      paste("bowtie -q -v 2 -m 1 --best --strata --sam -p", cores,
            spikeIdx,
            fq_unaligned,
            "| samtools view -bS -o ", bamSpike)
    }
  }, .(bamSpike, spikein_genome, fq_unaligned)]
  
  # UMI counts ----
  meta[, counts_cmd:= {
    if(overwrite | !file.exists(count))
      paste(Rpath,
            system.file("PROseq_pipeline", "count_UMIs.R", package = "vlfunctions"),
            bam,
            count)
  }, .(bam, count)]
  
  # SpikeIn counts ----
  meta[, spike_counts_cmd:= {
    if(overwrite | !file.exists(countSpike))
      paste(Rpath,
            system.file("PROseq_pipeline", "count_UMIs.R", package = "vlfunctions"),
            bamSpike,
            countSpike)
  }, .(bamSpike, countSpike)]
  
  # TSS count tables and stats ----
  meta[, count_tables_cmd:= {
    # Please specify:
    # [required] 1/ Reference genome count file
    # [required] 2/ Output folder (in which promoter/, gene_body/ and transcript/ sub-folders will be created)
    # [required] 3/ .rds file containing promoter annotations
    # [required] 4/ .rds file containing gene body annotations
    # [required] 5/ .rds file containing transcript annotations
    if(overwrite | any(!file.exists(c(countTables_promoter,
                                      countTables_geneBody,
                                      countTables_transcript))))
    {
      prom <- switch(genome,
                     "mm10"= "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_promoters.rds")
      geneBody <- switch(genome,
                         "mm10"= "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_genebody.rds")
      transcript <- switch(genome,
                           "mm10"= "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_transcript.rds")
      paste(Rpath,
            system.file("PROseq_pipeline", "compute_counts_table.R", package = "vlfunctions"),
            count, 
            paste0("db/count_tables/PROseq/", experiment, "/"),
            prom, 
            geneBody,
            transcript)
    }
  }, .(genome, experiment, count, countTables_promoter, countTables_geneBody, countTables_transcript)]
  
  # Read stats reference genome ----
  meta[, read_stats_cmd:= {
    # Please specify:
    # [required] 1/ Count file
    # [required] 2/ Output file
    if(overwrite | !file.exists(read_stats))
    {
      paste(Rpath,
            system.file("PROseq_pipeline", "compute_read_statistics.R", package = "vlfunctions"),
            count,
            read_stats)
    }
  }, .(count, read_stats)]
  
  # Read stats spike in ----
  meta[, spike_stats_cmd:= {
    # Please specify:
    # [required] 1/ Count file
    # [required] 2/ Output file
    if(overwrite | !file.exists(spike_stats))
    {
      paste(Rpath,
            system.file("PROseq_pipeline", "compute_read_statistics.R", package = "vlfunctions"),
            countSpike,
            spike_stats)
    }
  }, .(countSpike, spike_stats)]
  
  # Generate bw files ----
  meta[, bw_cmd:= {
    # Please specify:
    # [required] 1/ UMI counts file
    # [required] 2/ Output prefix (.ps.bw; .ns.bw)
    if(overwrite | any(!file.exists(c(bwPS, bwNS))))
    {
      paste(Rpath,
            system.file("PROseq_pipeline", "generate_bw_files.R", package = "vlfunctions"),
            count, 
            gsub(".ps.bw$", "", bwPS))
    }
  }, .(count, bwPS, bwNS)]
  
  # Return commands ----
  load_cmd <- paste(c(paste("cd", wdir),
                      "module load build-env/2020",
                      "module load cutadapt/1.18-foss-2018b-python-2.7.15",
                      "module load samtools/1.9-foss-2018b",
                      "module load bowtie/1.2.2-foss-2018b"), collapse = "; ")
  cols <- intersect(c("extract_cmd", "trim_cmd", "align_cmd", "unaligned_cmd", "alignSpike_cmd",
                      "counts_cmd", "spike_counts_cmd", "count_tables_cmd", "read_stats_cmd", "spike_stats_cmd", "bw_cmd"),
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
                name = "PROseq", 
                t = time,
                o= paste0(normalizePath(logs), "/", sampleID),
                e= paste0(normalizePath(logs), "/", sampleID))
      }, .(sampleID, cmd)] else
        return(cmd)
  }else
    warning("All output files already existed! No command submitted ;). Consider overwrite= T if convenient.")
}

#' DESeq2 analysis PROseq pipeline
#'
#' @param processed_metadata Path to the metadata file generated by vl_PROseq_processing, or the corresponding data.table.
#' @param FC_metadata_output An .rds path where to save the metadata file, which contains the directories containing .dds (DESeq2 objects) and FC table files and will be used to manage them.
#' By default, when the processed_metadata is a path to a processed_metadata files, "_FC_tables.rds" will be appended to the excel file path. 
#' @param Rpath Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".
#' @param dds_output_folder Output folder for .dds files (DESeq2 objects). Default= "db/dds/PROseq/".
#' @param FC_output_folder Output folder for FC tables. Default= "db/FC_tables/PROseq/".
#' @param PDF_output_folder Output folder for .pdf files of MA plots and statistics. Default= "pdf/PROseq/".
#' @param cores Number of cores per job. Default= 4.
#' @param mem Memory per job (in Go). Default= 16.
#' @param submit Should the command be submitted? default= FALSE.
#' @param wdir The working directory to use. Default= getwd(), meaning current working directory will be used.
#' @param logs Output folder for log files. Default= "db/logs/CUTNRUN/peak_calling".
#' @param time The time required for the SLURM scheduler. Default= '1-00:00:00'.
#'
#' @return Differential analyses using DESeq2.
#' @export
#'
#' @examples
#' # Make sure that all processed files exist:
#' processed <- readRDS("Rdata/metadata_PROseq_processed.rds")
#' file.exists(na.omit(unlist(processed[, fq1:bwNS])))
#' 
#' # Differential analysis (of note, no "overwrite" option is available in this case)
#' vl_PROseq_DESeq2(processed, # Here, you could also provide the path to the processed metadata table ("Rdata/metadata_PROseq_processed.rds").
#'                  FC_metadata_output = "Rdata/metadata_PROseq_FC_tables.rds",
#'                  submit = TRUE)
#'                  
#' # Once the commands ran, you can retrieve output folder using
#' FC <- readRDS("Rdata/metadata_PROseq_FC_tables.rds")
#' 
#' # For example, to retrieve spikein normalized FC tables:
#' files <- FC[, list.files(fcDir, paste0(DESeq2_condition, ".*", "spike")), .(fcDir, DESeq2_condition)]$V1
#'                  
vl_PROseq_DESeq2 <- function(processed_metadata, ...) UseMethod("vl_PROseq_DESeq2")

#' @describeIn vl_PROseq_DESeq2  for processed_metadata file path
#' @export
vl_PROseq_DESeq2.character <- function(processed_metadata,
                                       FC_metadata_output= gsub("processed.rds$", "FC_tables.rds", processed_metadata),
                                       ...)
{
  meta <- readRDS(processed_metadata)
  vl_PROseq_DESeq2(processed_metadata= meta,
                   FC_metadata_output= FC_metadata_output,
                   ...)
}

#' @describeIn vl_PROseq_DESeq2 default method
#' @export
vl_PROseq_DESeq2.default <- function(processed_metadata,
                                     FC_metadata_output,
                                     Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                     dds_output_folder= "db/dds/PROseq/",
                                     FC_output_folder= "db/FC_tables/PROseq/",
                                     PDF_output_folder= "pdf/PROseq/",
                                     cores= 4,
                                     mem= 16,
                                     submit= FALSE,
                                     wdir= getwd(),
                                     logs= "db/logs/PROseq/DESeq2",
                                     time= '1-00:00:00')
{
  # Hard copy metadata ----
  meta <- data.table::copy(processed_metadata)
  
  # Melt features ----
  meta <- melt(meta,
               id.vars = c("sampleID", "experiment", "DESeq2_name", "DESeq2_condition", "DESeq2_control", "read_stats", "spike_stats"),
               measure.vars = patterns("countTable"),
               value.name = "countFile",
               variable.name = "feature")
  meta[, feature:= gsub(".*_(.*$)", "\\1", feature)]
  meta[feature=="geneBody", feature:= "gene_body"]
  
  # Create output directories/prefixes ----
  meta[, ddsPrefix:= paste0(dds_output_folder, experiment, "/", feature)]
  meta[, ddsStats:= paste0(PDF_output_folder, experiment, "_", feature, "_reads_statistics.pdf")]
  meta[, fcDir:= paste0(FC_output_folder, experiment, "/", feature, "/")]
  meta[, MAplot_prefix:= paste0(PDF_output_folder, "MA_plots/", experiment, "_", feature)]
  
  # Save processed FC metadata ----
  if(!grepl(".rds$", FC_metadata_output))
    stop("FC_metadata_output should end up with a .rds extension") else
      saveRDS(meta, FC_metadata_output)
  
  # Create directories ----
  dirs <- dirname(unlist(meta[, .(ddsPrefix, ddsStats, MAplot_prefix)]))
  dirs <- na.omit(unique(c(meta$fcDir, dirs, logs)))
  if(any(!dir.exists(dirs)))
  {
    sapply(dirs, dir.create, showWarnings = F, recursive = T)
    outDirs <- paste0(c(FC_output_folder, PDF_output_folder), collapse= ", ")
    print(paste("Output directories were created in:", outDirs, "!"))
  }
  
  # DESeq2 commands ---- 
  meta[, DESeq2_cmd:= {
    # Each command produces a unique dds file and all the related FC tables/MA plots
    # [required] 1/ A comma-separated list of count (ref genome)
    # [required] 2/ A comma-separated list of read statistics (reference genome) 
    # [required] 3/ A comma-separated list of spike-in statistics 
    # [required] 4/ A comma-separated list of sample names 
    # [required] 5/ A comma-separated list of conditions 
    # [required] 6/ A comma-separated list of controls 
    # [required] 7/ .dds output prefix 
    # [required] 8/ dds statistics output 
    # [required] 9/ FC tables output folder 
    # [required] 10/ MAplot prefix
    paste(Rpath,
          system.file("PROseq_pipeline", "DESeq2_analysis.R", package = "vlfunctions"),
          paste0(countFile, collapse = ","), 
          paste0(read_stats, collapse = ","),
          paste0(spike_stats, collapse = ","),
          paste0(DESeq2_name, collapse = ","),
          paste0(DESeq2_condition, collapse = ","),
          paste0(DESeq2_control, collapse = ","),
          ddsPrefix,
          ddsStats,
          fcDir,
          MAplot_prefix)
  }, .(experiment, ddsPrefix, ddsStats, fcDir, MAplot_prefix)]
  
  # Submit commands ----
  cmd <- meta[, .(cmd= paste0(unique(na.omit(unlist(.SD))), collapse= "; ")), .(experiment, feature), .SDcols= patterns("_cmd$")]
  if(nrow(cmd))
  {
    if(submit)
      cmd[, {
        vl_bsub(cmd, 
                cores= cores, 
                m = mem, 
                name = "PROseq", 
                t = time,
                o= paste0(normalizePath(logs), "/", paste0(experiment, "_", feature)),
                e= paste0(normalizePath(logs), "/", paste0(experiment, "_", feature)))
      }, .(experiment, feature, cmd)] else
        return(cmd)
  }else
    warning("No commands detected. Is the metadata empty?")
}