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

#' Submits to R using singularity
#' @export
vl_Rsub_singularity <- function(R_script, args_v= NULL)
{
  args_v <- paste0(args_v, collapse= " ")
  cmd <- paste0("module load singularity/3.4.1;
         /opt/ohpc/pub/libs/singularity/3.4.1/bin/singularity run --app Rscript /groups/stark/software-all/singularity_img/singularity.R.with_pckg.simg ", R_script)
  if(!is.null(args_v))
    cmd <- paste0(cmd, " ", args_v)
  return(cmd)
}

#' Wrap commands to send to bash on windows
#' @export

vl_bash_wrap_windows <- function(cmd)
{
  cmd <- paste("bash -c \"", cmd, "\"")  
}

#' Extract reads from bam VBC
#' @export
vl_extract_reads_VBC <- function(bam, output_prefix, rev_comp_i5)
{
  paste0("module load build-env/2020; module load samtools/1.9-foss-2018b; 
         /groups/stark/software-all/shell/demultiplexPE.sh -i ", bam, " -o ", output_prefix, " -b ", '"', rev_comp_i5, '"', " -u TRUE")
}

#' BOWTIE Alignment to custom genome
#' @export
vl_bowtie1_align <- function(fq1, index_prefix, sam_output, o= NA, e= NA, execute= F, cores= 6, d= NA)
{
  paste0("module load bowtie/1.2.2-foss-2018b; /software/2020/software/bowtie/1.2.2-foss-2018b/bin/bowtie ", 
         index_prefix, " -p ", cores, " -m1 ", fq1, " -v 2 -m 1 --best --strata -S ", sam_output)
}

#' BOWTIE2 Alignment to custom genome
#' @export
vl_bowtie2_align <- function(fq1, fq2= NA, index_prefix, sam_output, cores= 6)
{
  if(is.na(fq2))
    files <- paste0(" -U ", fq1) else
      files <- paste0(" -1 ", fq1, " -2 ", fq2)
  paste0("module load bowtie2/2.3.4.2-foss-2018b; /software/2020/software/bowtie2/2.3.4.2-foss-2018b/bin/bowtie2 -p ", 
         cores, " -x ", index_prefix, files, " -S ", sam_output)
}

# Trim 3' fastq reads using trimgalore
#' @export
vl_trim3_fastq <- function(fq, keep, output_folder, cores= 1)
{
  paste0("module load trim_galore/0.6.2-foss-2018b-python-3.6.6; /software/2020/software/trim_galore/0.6.2-foss-2018b-python-3.6.6/trim_galore --gzip --hardtrim3 ", keep, " --cores ", cores, " -o ", output_folder, " ", fq)
}

# Trim 5' fastq reads using trimgalore
#' @export
vl_trim5_fastq <- function(fq, keep, output_folder, cores= 1)
{
  paste0("module load trim_galore/0.6.2-foss-2018b-python-3.6.6; /software/2020/software/trim_galore/0.6.2-foss-2018b-python-3.6.6/trim_galore --gzip --hardtrim5 ", keep, " --cores ", cores, " -o ", output_folder, " ", fq)
}

# bamtobed wrapper
#' @export
vl_bamToBed <- function(bam, bed_output, cigar= F)
{
  cmd <- "/software/2020/software/bedtools/2.27.1-foss-2018b/bin/bamToBed"
  if(cigar)
    cmd <- paste0(cmd, " -cigar")
  paste0(cmd, " -i ", bam, " > ", bed_output)
}

# bamtobedpe wrapper
#' @export
vl_bamToBedpe <- function(bam, bed_output, cigar= F)
{
  cmd <- "/software/2020/software/bedtools/2.27.1-foss-2018b/bin/bamToBed"
  if(cigar)
    cmd <- paste0(cmd, " -cigar")
  paste0(cmd, " -bedpe -i ", bam, " > ", bed_output)
}

# intersectbed wrapper
#' @export
vl_intersectBed <- function(a, b, same_strand=F, sorted= F)
{
  cmd <- paste0("/software/2020/software/bedtools/2.27.1-foss-2018b/bin/intersectBed -c ")
  if(sorted)
    cmd <- paste0(cmd, "-sorted ")  
  if(same_strand)
    cmd <- paste0(cmd, "-s ")
  paste0(cmd, "-a ", a, " -b ", b)
}

# closestbed wrapper
#' @export
vl_closestBed <- function(a, b, k= 1)
{
  paste0("/software/2020/software/bedtools/2.27.1-foss-2018b/bin/closestBed -d -k ", k, " -a ", a, " -b ", b)
}

# bigBedToBed wrapper
#' @export
vl_bigBedToBed <- function(bb, bed)
{
  paste0("/software/2020/software/kent_tools/20190507-linux.x86_64/bin/bigBedToBed ", bb, " ", bed)
}

# samtoBam wrapper
#' @export
vl_samtools_samToBam <- function(sam_bam, bam_output, cores= 6)
{
  paste0("module load samtools/1.9-foss-2018b;
         /software/2020/software/samtools/1.9-foss-2018b/bin/samtools view -@ ", cores-1, " -b -o ", bam_output, " ", sam_bam)
}

# mapq filter sam wrapper
#' @export
vl_samtools_mapq_filter <- function(sam_bam, bam_output, q= 30, cores= 6)
{
  paste0("module load samtools/1.9-foss-2018b;
         /software/2020/software/samtools/1.9-foss-2018b/bin/samtools view -@ ", cores-1, " -b -q ", q, " -o ", bam_output, " ", sam_bam)
}

# Sort bam file
#' @export
vl_samtools_bam_sort <- function(sam_bam, bam_output, cores= 6)
{
  # Sort bam
  paste0("module load samtools/1.9-foss-2018b;
         /software/2020/software/samtools/1.9-foss-2018b/bin/samtools sort -n -@ ", cores-1, " -o ", bam_output, " ", sam_bam)
}

# Check if bam file sorted
#' @export
vl_samtools_check_if_sorted <- function(sam_bam)
{
  # Check bam
  paste0("module load samtools/1.9-foss-2018b;
         /software/2020/software/samtools/1.9-foss-2018b/bin/samtools view -H ", sam_bam, " | grep @HD")
}
