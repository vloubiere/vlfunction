#' Stark lab wrappers
#'
#' bsub submits to the cluster
#' 
#' @param cmd Command to be submitted
#' @param cores Number of cores. default= 12L
#' @param m Memory in Go. Default= 12
#' @param t Time parameter for slurm in the format hh:mm:ss. Default= '08:00:00'. For 2 days:'2-00:00:00'
#' @param name Name of the job. Default= "vl
#' @param o Std output directory
#' @param e Std error directory
#' @param wdir Working directory. defaults to getwd()
#' @param execute If FALSE, simply returns the full command. Default= TRUE
#' @param d Name of te dependent job
#'
#' @export
vl_bsub <- function(cmd,
                    cores= 6L,
                    m= 12,
                    t= '08:00:00',
                    name= "vl",
                    o= NA,
                    e= NA,
                    wdir= getwd(),
                    execute= T,
                    d= NA)
{
  # DO NOT MODIFY UNLESS YOU EXACTLY KNOW WHAT YOU ARE DOING (all other functions use it)
  # gen cmd
  bsub_cmd <- paste0("/groups/stark/software-all/shell/bsub_gridengine -C ", cores , " -n ", name)
  if(any(!is.na(d)))
    bsub_cmd <- paste0(bsub_cmd, " -d ", d[1])
  if(length(d)>1)
    bsub_cmd <- paste0(bsub_cmd, paste0(",afterok:", d[-1], collapse= ""))
  if(!is.na(m)){bsub_cmd <- paste0(bsub_cmd, " -m ", m)}
  if(!is.na(o)){bsub_cmd <- paste0(bsub_cmd, " -o ", o)}
  if(!is.na(e)){bsub_cmd <- paste0(bsub_cmd, " -e ", e)}
  if(!is.na(t)){bsub_cmd <- paste0(bsub_cmd, " -T ", t)}
  if(!is.na(t)){bsub_cmd <- paste0(bsub_cmd, " -T ", t)}
  # cd to wdir
  if(!is.na(wdir)){cmd <- paste0("cd ", wdir, "; ", cmd)}
  bsub_cmd <- paste0(bsub_cmd, " \"", cmd, "\"")
  
  # Submit
  if(isTRUE(execute))
  {
    Sys.unsetenv("SBATCH_RESERVATION")
    Sys.unsetenv("SBATCH_WCKEY")
    # Write command in file
    tmp <- tempfile(tmpdir = "/scratch/stark/vloubiere/",
                    fileext = ".sh")
    writeLines(bsub_cmd, tmp)
    # Execute file using ssh
    job_ID <- system(paste0("ssh localhost sh ", tmp),
                     intern = T,
                     ignore.stderr = TRUE)
    file.remove(tmp)
    return(gsub(".* (.*$)", "\\1", job_ID[2]))
  }else
  {
    return(bsub_cmd)
  }
}

#' Submits to R on cluster
#' @export
vl_Rsub <- function(R_script, args_v= NULL)
{
  # Submit R job wrapper ####
  args_v <- paste0(args_v, collapse= " ")
  cmd <- paste0("module load r/3.6.2-foss-2018b; /software/2020/software/r/3.6.2-foss-2018b/bin/Rscript ", R_script)
  if(!is.null(args_v))
    cmd <- paste0(cmd, " ", args_v)
  return(cmd)
}

#' Extract reads from bam VBC
#' @export
vl_extract_reads_VBC <- function(bam,
                                 output_prefix,
                                 rev_comp_i5)
{
  paste0("module load build-env/2020; module load samtools/1.9-foss-2018b; 
         /groups/stark/software-all/shell/demultiplexPE.sh -i ", bam, " -o ", output_prefix, " -b ", '"', rev_comp_i5, '"', " -u TRUE")
}

#' Use dropbox uploader
#' @export
vl_dropbox_upload <- function(local_path,
                              remote_path)
{
  cmd <- paste("sh  /groups/stark/vloubiere/apps/dropbox/dropbox_uploader.sh upload", local_path, remote_path)
  system(cmd)
}

#' Use dropbox downloader
#' @export
vl_dropbox_download <- function(local_path,
                                remote_path)
{
  cmd <- paste("sh  /groups/stark/vloubiere/apps/dropbox/dropbox_uploader.sh download", remote_path, local_path)
  system(cmd)
}

#' Download iCisTarget report
#' 
#' Download an iCisTarget report, uncompresses it and stores it in a subfolder
#'
#' @param url URL where the reports is to be found
#' @param dir Directory where the report should be saved. 
#' @param folder_name Name of the unziped subfolder. By default, the downloaded subfolder will be named using the iCisTarget job name.
#'
#' @export
vl_download_iCistarget <- function(url, 
                                   dir,
                                   folder_name)
{
  cmd <- paste0("cd ", normalizePath(dir), ";wget ", url, "; unzip archive.zip; rm archive.zip")
  system(cmd)
  if(missing(folder_name))
    folder_name <- unlist(fread(paste0(normalizePath(dir), "/icistarget/statistics.tbl"))[1, 1])
  final <- paste0(normalizePath(dir), "/", folder_name)
  cmd <- paste("mv", paste0(normalizePath(dir), "/icistarget"), final)
  system(cmd)
}

# Squeue to the system
#' @export
vl_squeue <- function()
{
  system("ssh localhost squeue -u vincent.loubiere",
         ignore.stderr = T)
}

# Scancel all jobs except Rstudio
#' @export
vl_scancel <- function()
{
  .c <- fread(cmd= "ssh localhost squeue -u vincent.loubiere")
  .c <- .c[!(NAME %in% c("[RStudio", "jupyter_"))]
  if(nrow(.c))
    system(paste(c("ssh localhost scancel", .c$JOBID), collapse= " "),
           ignore.stderr = T)
}

#' Import bam file with all fields
#'
#' @param bam_path 
vl_import_bam <- function(bam_path,
                          ncores= data.table::getDTthreads()-1)
{
  cmd <- "module load build-env/2020; module load samtools/1.9-foss-2018b; module load bedtools/2.17.0-foss-2018b; samtools view -@"
  cmd <- paste(cmd, ncores, bam_path)
  bam <- fread(cmd= cmd, fill = T)
  return(bam)
}

#' Submits to R using singularity
vl_Rsub_singularity <- function(R_script,
                                args_v= NULL)
{
  args_v <- paste0(args_v, collapse= " ")
  cmd <- paste0("module load singularity/3.4.1;
         /opt/ohpc/pub/libs/singularity/3.4.1/bin/singularity run --app Rscript /groups/stark/software-all/singularity_img/singularity.R.with_pckg.simg ", R_script)
  if(!is.null(args_v))
    cmd <- paste0(cmd, " ", args_v)
  return(cmd)
}

#' See most recent error file in a directory
#' 
#' @param dir URL where the reports is to be found
#'
#' @export
vl_last_err <- function(dir)
{
  dat <- data.table(file= list.files(dir, ".err$", full.names = T))
  dat[, mtime:= file.info(file)$mtime]
  last <- dat[which.max(mtime)]
  print(paste("Printed on:", last$mtime))
  file.show(last$file)
}

#' See most recent error file in a directory
#' 
#' @param metadata Path to an excel metadata sheet
#' @param first.col.name If metadata in in .xlsx format, this character string will be used to detect the starting line. Default= "user".
#'
#' @export
vl_import_metadata_sheet <- function(metadata, first.col.name= "user")
{
  meta <- if(grepl(".txt$", metadata))
  {
    return(fread(metadata))
  }else if(grepl(".rds$", metadata))
  {
    return(readRDS(metadata))
  }else if(grepl(".xlsx$", metadata))
  {
    sheet <- readxl::read_xlsx(metadata)
    # Slip lines before first column names
    sel.rows <- cumsum(c(sheet[[1]]==first.col.name))>0
    if(sum(sel.rows, na.rm= TRUE)==0)
      stop("The 'user' column, which should be the first of the excel sheet, could not be found!")
    sheet <- sheet[(sel.rows),]
    meta <- as.data.table(sheet[-1,])
    names(meta) <- unlist(sheet[1,])
    return(meta)
  }else
    stop("metadata file extension should be one of '.txt', '.rds' or '.xlsx'!")
}
