#' Stark lab wrappers
#'
#' bsub submits to the cluster
#' @export

vl_bsub <- function(cmd, cores= 1, name= "vl", d= NA, m= NA, o= NA, e= NA, t= NA, execute= T)
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
  bsub_cmd <- paste0(bsub_cmd, " \"", cmd, "\"")
  
  # Submit
  if(isTRUE(execute))
  {
    Sys.unsetenv("SBATCH_RESERVATION")
    Sys.unsetenv("SBATCH_WCKEY")
    job_ID <- system(bsub_cmd, intern= T)
    return(unlist(data.table::tstrsplit(job_ID[2], " ", keep= 4)))
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
vl_extract_reads_VBC <- function(bam, output_prefix, rev_comp_i5)
{
  paste0("module load build-env/2020; module load samtools/1.9-foss-2018b; 
         /groups/stark/software-all/shell/demultiplexPE.sh -i ", bam, " -o ", output_prefix, " -b ", '"', rev_comp_i5, '"', " -u TRUE")
}

#' Use dropbox uploader
#' @export
vl_dropbox_upload <- function(local_path, remote_path)
{
  cmd <- paste("sh  /groups/stark/vloubiere/apps/dropbox/dropbox_uploader.sh upload", local_path, remote_path)
  system(cmd)
}

#' Use dropbox downloader
#' @export
vl_dropbox_download <- function(local_path, remote_path)
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

#' Import paired bam file as bedpe
#'
#' @param bam_path 
#' @param bedpe Imports bam as BEDPE format. Requires BAM to be grouped or sorted by query.
#' @param ed Use BAM edit distance (NM tag) for BED score. Default for BED is to use mapping quality. Default for BEDPE is to use the minimum of the two mapping qualities for the pair. When -ed is used with bedpe= T, the total edit distance from the two mates is reported.
vl_import_bamToBed <- function(bam_path, bedpe= F, ed= F)
{
  cmd <- "module load build-env/2020; module load bedtools/2.17.0-foss-2018b; bamToBed -i"
  cmd <- paste(cmd, bam_path)
  if(bedpe)
    cmd <- paste(cmd, "-bedpe")
  if(ed)
    cmd <- paste(cmd, "-ed")
  bed <- fread(cmd= cmd)
  return(bed)
}

#' Import bam file with all fields
#'
#' @param bam_path 
vl_import_bam <- function(bam_path, ncores= data.table::getDTthreads()-1)
{
  cmd <- "module load build-env/2020; module load samtools/1.9-foss-2018b; module load bedtools/2.17.0-foss-2018b; samtools view -@"
  cmd <- paste(cmd, ncores, bam_path)
  bam <- fread(cmd= cmd, fill = T)
  return(bam)
}

#' Submits to R using singularity
vl_Rsub_singularity <- function(R_script, args_v= NULL)
{
  args_v <- paste0(args_v, collapse= " ")
  cmd <- paste0("module load singularity/3.4.1;
         /opt/ohpc/pub/libs/singularity/3.4.1/bin/singularity run --app Rscript /groups/stark/software-all/singularity_img/singularity.R.with_pckg.simg ", R_script)
  if(!is.null(args_v))
    cmd <- paste0(cmd, " ", args_v)
  return(cmd)
}