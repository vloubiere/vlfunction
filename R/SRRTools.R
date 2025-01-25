#' Title
#'
#' @param SRR 
#' @param layout 
#' @param output_folder 
#' @param tmp_folder 
#' @param threads 
#' @param mem 
#'
#' @return
#' @export
#'
#' @examples
vl_download_SRA <- function(SRR,
                            layout,
                            output_folder,
                            tmp_folder,
                            threads= 8,
                            mem= 2048)
{
  if(!dir.exists(output_folder))
    dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(tmp_folder))
    dir.create(tmp_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Initiate command ----
  cmd <- "module load build-env/f2022; module load sra-toolkit/3.1.1-centos_linux64; fasterq-dump"
  if(layout=="PAIRED")
    cmd <- paste(cmd, "--split-files")
  cmd <- paste(cmd, "--threads", threads, "--mem", mem, "--temp", tmp_folder, "--outdir", output_folder, SRR)
  
  # Return
  return(cmd)
}