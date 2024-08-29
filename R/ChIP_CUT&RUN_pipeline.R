#' ORFtag pipeline
#' 
#' Takes as input a (corerectly formated) metadata file, save the processed metadata file and returns the command lines to 1/ extract reads from VBC bam file, 2/ trim the reads and align to mouse/human genome (see 'species' column of the metadata table) and return alignment statistics as well as collapsed reads and 4/ assign insertions to closest downstream genes.
#'
#' @param metadata The path to a .txt, tab-separated metadata file containing at least 12 columns. See vlfunctions::vl_metadata_ORFtag for an example.
#' @param processed_metadata_output An .rds path where to save the processed metadata file (containing the paths of output files). By default, "_processed.rds" will be appended to the metadata path.
#' @param scratch_folder Folder where intermediate bam and fastq files will be saved
#' @param cores Number of cores per job. Default= 8
#' @param mem Memory per job (in Go). Default= 32.
#' @param overwrite Should existing files be overwritten?
#' @param submit Should the command be submitted? default= G
#' @param wdir The working directory to use. defaut= getwd().
#' @param logs Path to save logs. Default= "db/logs/CUTNRUN/processing"
#' @param time The time required for the SLURM scheduler. Default= '1-00:00:00'
#'
#' @return Return a data.table containing, for each sampleID, the concatenated commands that are required to process the data. These commands can then be submitted using ?vl_bsub().
#' @export
#'
#' @examples
#' Takes as input a .xlsx metadata sheet that should be formatted as follows:
#' >>>>
#' READ ME																									
#' sampleID should be the concatenation of ab, cell_line, treatment, experiment and replicate; and should be unique for each sample/replicate. Samples with the same sampleID will be merged (re-sequencing...).																									
#' The condition should be shared between replicates, and will be used to call peaks on the merge data.																									
#' The input column should specify the replicate-matched input sampleID.																									
#' user	method	experiment	antibody	cell_line	treatment	replicate	condition	sampleID	barcode	i5	genome	input	layout	sequencer	bam_path
#' <<<<
#' custom extra columns are tolerated						
#' 
#' Print command lines for all samples without submitting the jobs:
#' check <- vl_CUTNRUN_processing(metadata = "/groups/stark/nemcko/Hcfc1_paper/Hcfc1_new/Rdata/metadata_CutNRun.xlsx",
#'                                processed_metadata_output = "Rdata/CUTNRUN_metadata_processed.rds",
#'                                overwrite = T,
#'                                submit = F)
#' 
#' Run the pipeline without overwriting existing files:
#' vl_CUTNRUN_processing(metadata = "/groups/stark/nemcko/Hcfc1_paper/Hcfc1_new/Rdata/metadata_CutNRun.xlsx",
#'                       processed_metadata_output = "Rdata/CUTNRUN_metadata_processed.rds",
#'                       overwrite = F,
#'                       submit = T)
#'                       
vl_CUTNRUN_processing <- function(metadata,
                                  processed_metadata_output= gsub(".txt$", "_processed.rds", metadata),
                                  scratch_folder= "/scratch/stark/vloubiere/CUTNRUN",
                                  cores= 8,
                                  mem= 64,
                                  overwrite= F,
                                  submit= F,
                                  wdir= getwd(),
                                  logs= "db/logs/CUTNRUN/processing",
                                  time= '1-00:00:00')
{
  # Import metadata and check format ----
  meta <- readxl::read_xlsx(metadata, skip = 4)
  meta <- as.data.table(meta)
  cols <- c("user","method","experiment","antibody","cell_line","treatment","replicate","condition","sampleID","barcode","i5","genome","input","layout","sequencer","bam_path")
  if(!all(cols %in% names(meta)))
    stop(paste(c("metadata file should contain at least these 12 columns:", cols), collapse = " "))
  if(any(!na.omit(meta$input) %in% meta$sampleID))
    stop("Some input(s) have no correspondence in the sampleID(s) column!")
  if(any(!meta$genome %in% c("mm10", "hg38")))
    stop("For now, only supported genomes are mm10 and hg38!")
  if(!all(meta[, sampleID==paste0(antibody, "_", cell_line, "_", treatment, "_", experiment, "_", replicate)]))
    warning("sampleID should be the catenation of Ab, Cell_line, Treatment, experiment and replicate. Any samples with the same sampleID will be collapsed!")
  
  # Generate output paths ----
  meta[(layout=="SINGLE"), fq1:= paste0(scratch_folder, "/fq/", gsub(".bam", "", basename(bam_path)), "_", make.unique(sampleID), ".fq.gz")]
  meta[(layout=="PAIRED"), fq1:= paste0(scratch_folder, "/fq/", gsub(".bam", "", basename(bam_path)), "_", make.unique(sampleID), "_1.fq.gz")]
  meta[(layout=="PAIRED"), fq2:= gsub("_1.fq.gz$", "_2.fq.gz", fq1)]
  meta[, fq1_trim:= gsub(".fq.gz$", ifelse(layout=="PAIRED", "_val_1.fq.gz", "_trimmed.fq.gz"), fq1), layout]
  meta[, fq2_trim:= gsub(".fq.gz$", "_val_2.fq.gz", fq2)]
  meta[, bam:= paste0(scratch_folder, "/bam/", sampleID, ".bam")] # re-sequencing are merged from this step on!
  meta[meta, bam_input:= i.bam, on= "input==sampleID"] # Input bam
  meta[, alignment_stats:= paste0("db/alignment_stats/CUTNRUN/", gsub(".bam$", "_stats.txt", basename(bam)))]
  meta[, mapq30_stats:= gsub("_stats.txt$", "_mapq30_stats.txt", alignment_stats)]
  meta[, peaks_file:= paste0("db/peaks/CUTNRUN/replicates/", sampleID, "_peaks.narrowPeak")]
  meta[, bw_file:= paste0("db/bw/CUTNRUN/replicates/", sampleID, ".bw")]
  meta[, twoReps:= length(unique(bam))>1, condition] # Check if at least two replicates for conf. peaks
  meta[(twoReps), merged_peaks_file:= paste0("db/peaks/CUTNRUN/merge/", condition, "_merged_peaks.narrowPeak")]
  meta[(twoReps), merged_bw_file:= paste0("db/bw/CUTNRUN/merge/", condition, "_merged.bw")]
  meta[(twoReps), confident_peaks_file:= paste0("db/peaks/CUTNRUN/confident/", condition, "_confident_peaks.narrowPeak")]
  
  # Save processed metadata ----
  if(!grepl(".rds$", processed_metadata_output))
    stop("processed_metadata_output should end up with a .rds extension") else
      saveRDS(meta, processed_metadata_output)
  
  # Create output directories ----
  dirs <- c(logs,
            na.omit(unique(dirname(unlist(meta[, fq1:confident_peaks_file])))))
  if(any(!dir.exists(dirs)))
  {
    sapply(dirs, dir.create, showWarnings = F, recursive = T)
    print("Output directories were created in db/ and scratch/!")
  }
  
  # Print conditions ----
  cat(paste(length(unique(meta$condition)), "conditions detected, of which", length(unique(meta[(twoReps), condition])), "had >1 replicates:\n"))
  meta[(twoReps), cat(paste0(condition, ":\n\t", paste0(unique(sampleID), collapse= "\n\t"), "\n")), condition]
  cat(paste(length(unique(meta[(!twoReps), condition])), "samples had only one replicate:"))
  meta[(!twoReps), cat(paste0(condition, ":\n\t", paste0(unique(sampleID), collapse= "\n\t"), "\n")), condition]
  
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
                   "; gzip", gsub(".gz$", "", fq1))
      if(!is.na(fq2))
        cmd <- paste0(cmd, "; gzip ", gsub(".gz$", "", fq2))
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
            "| samtools sort -@", cores-1, "-T", paste0(scratch_folder, "/bam/sorted"), # Sort with temp files in scratch folder
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
#' @param processed_metadata Processed metadata file generated by vl_CUTNRUN_processing
#' @param Rpath The path to the Rscript executable to use, on which the latest version of the vlfunction package should be installed. Default= /software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript
#' @param extsise The extsize to be used, meaning that the building of the peak model will be skipped (--nomodel). To build the model, set extize= NA. Default= 300 
#' @param cores Number of cores per job. Default= 8
#' @param mem Memory per job (in Go). Default= 32.
#' @param overwrite Should existing files be overwritten?
#' @param submit Should the command be submitted? default= G
#' @param wdir The working directory to use. defaut= getwd().
#' @param logs Path to save logs. Default= "db/logs/CUTNRUN/peak_calling"
#' @param time The time required for the SLURM scheduler. Default= '1-00:00:00'
#'
#' @return Command lines for peak calling.
#' @export
#'
#' @examples
#' Takes as input a processed metadata file produced by ?vl_CUTNRUN_processing().
#' 
#' Print the commands without submitting:
#' check <- vl_CUTNRUN_peakCalling(processed_metadata = "Rdata/CUTNRUN_metadata_processed.rds",
#'                                 overwrite = T,
#'                                 submit = F)
#' Run the commands without overwriting existing files:
#' vl_CUTNRUN_peakCalling(processed_metadata = "Rdata/CUTNRUN_metadata_processed.rds",
#'                        overwrite = T,
#'                        submit = T)
#' 
vl_CUTNRUN_peakCalling <- function(processed_metadata,
                                   Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                                   extsize= 300,
                                   cores= 8,
                                   mem= 64,
                                   overwrite= F,
                                   submit= F,
                                   wdir= getwd(),
                                   logs= "db/logs/CUTNRUN/peak_calling",
                                   time= '1-00:00:00')
{
  # Import metadata and check format ----
  meta <- readRDS(processed_metadata)[!is.na(input)]
  
  # Create missing output dirs ----
  if(!dir.exists(logs))
  {
    dir.create(logs, showWarnings = F, recursive = T)
    print(paste(logs, "directory was created!"))
  }
  
  # Check if at least two replicates ----
  meta[(!twoReps), print(paste("No replicates found for sample", sampleID, "-> confident peaks calling skipped")), sampleID]

  # Peak calling on individual replicates ----
  meta[, peak_cmd:= {
    if(overwrite | !file.exists(peaks_file))
    {
      gAbb <- fcase(genome=="mm10", "mm",
                    genome=="hg38", "hg")
      paste("macs2 callpeak -B --SPMR --keep-dup 1 -g", gAbb, 
            ifelse(layout=="PAIRED", "-f BAMPE", ""),
            ifelse(is.na(extsize), "", paste("--nomodel --extsize", extsize)),
            "--outdir", paste0(dirname(peaks_file), "/"),
            "-t", bam,
            "-c", bam_input,
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
      gAbb <- fcase(genome=="mm10", "mm",
                    genome=="hg38", "hg")
      paste("macs2 callpeak -B --SPMR --keep-dup 1 -g", gAbb,
            ifelse(layout=="PAIRED", "-f BAMPE", ""),
            ifelse(is.na(extsize), "", paste("--nomodel --extsize", extsize)),
            "--outdir", paste0(dirname(merged_peaks_file), "/"),
            "-t", paste0(unique(bam), collapse= " "),
            "-c", paste0(unique(bam_input), collapse= " "),
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
