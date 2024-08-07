% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ORFtag_pipeline.R
\name{vl_ORFtag_pipeline}
\alias{vl_ORFtag_pipeline}
\title{ORFtag pipeline}
\usage{
vl_ORFtag_pipeline(
  metadata,
  processed_metadata_output,
  Rpath = "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
  cores = 8,
  mem = 32,
  overwrite = F,
  submit = F,
  wdir = getwd(),
  logs = "db/logs/ORFtag"
)
}
\arguments{
\item{metadata}{The path to a .txt, tab-separated metadata file containing at least 12 columns. See vlfunctions::vl_metadata_ORFtag for an example.}

\item{processed_metadata_output}{An .rds path where to save the processed metadata file (containing the paths of output files).}

\item{Rpath}{Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript"}

\item{cores}{Number of cores per job. Default= 8}

\item{mem}{Memory per job (in Go). Default= 32.}

\item{overwrite}{Should existing files be overwritten?}

\item{submit}{Should the command be submitted? default= G}

\item{wdir}{The working directory to use. defaut= getwd().}

\item{logs}{Path to save logs. Default= "db/logs"}
}
\value{
Return a data.table containing, for each sampleID, the concatenated commands that are required to process the data. These commands can then be submitted using ?vl_bsub().
}
\description{
Takes as input a (corerectly formated) metadata file, save the processed metadata file and returns the command lines to 1/ extract reads from VBC bam file, 2/ trim the reads and align to mouse/human genome (see 'species' column of the metadata table) and return alignment statistics as well as collapsed reads and 4/ assign insertions to closest downstream genes.
}
\examples{
path <- "Rdata/metadata_ORFtag.txt"
fwrite(vlfunctions::vl_metadata_ORFtag, path)

# Process data ----
vl_ORFtag_pipeline(metadata= path,
                   processed_metadata_output = "Rdata/metadata_ORFtag_processed.rds",
                   Rpath= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
                   cores= 8,
                   mem= 32,
                   overwrite= F,
                   submit= T,
                   logs= "db/logs")

# Call hits ----
vl_ORFtrap_call_hits(sorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_sort_rep1_same_strand.txt",
                                               "db/gene_assignment/ORFtag/Activator2_sort_rep2_same_strand.txt"),
                     unsorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_input_rep1_same_strand.txt",
                                                 "db/gene_assignment/ORFtag/Activator2_input_rep2_same_strand.txt"),
                     exons_start_gtf = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
                     name = "Activator_2",
                     output_suffix = "_vs_input.txt",
                     output_folder_FC_file = "db/FC_tables/ORFtag/")

# Call hits using revese strand (sanity check -> be cautious with the hits that are also found here!) ----
vl_ORFtrap_call_hits(sorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_sort_rep1_rev_strand.txt",
                                               "db/gene_assignment/ORFtag/Activator2_sort_rep2_rev_strand.txt"),
                     unsorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_input_rep1_rev_strand.txt",
                                                 "db/gene_assignment/ORFtag/Activator2_input_rep2_rev_strand.txt"),
                     exons_start_gtf = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
                     name = "Activator_2",
                     output_suffix = "_vs_input_rev_strand.txt",
                     output_folder_FC_file = "db/FC_tables/ORFtag/")

hits <- fread("db/FC_tables/ORFtag/Activator_2_vs_input.txt")[(hit), gene_name]
sanityCheck <- fread("db/FC_tables/ORFtag/Activator_2_vs_input_rev_strand.txt")[(hit), gene_name]
print(paste(length(hits), "hits were found, out of which", sum(hits \%in\% sanityCheck),
            "were also found using reversed strand and should be considered carefully"))

# Call hits using strand bias (not used) ----
vl_ORFtrap_call_hits_strandBias(sorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_sort_rep1_same_strand.txt",
                                                          "db/gene_assignment/ORFtag/Activator2_sort_rep2_same_strand.txt"),
                                sorted_reverse_counts = c("db/gene_assignment/ORFtag/Activator2_sort_rep1_rev_strand.txt",
                                                          "db/gene_assignment/ORFtag/Activator2_sort_rep2_rev_strand.txt"),
                                unsorted_forward_counts = c("db/gene_assignment/ORFtag/Activator2_input_rep1_same_strand.txt",
                                                            "db/gene_assignment/ORFtag/Activator2_input_rep2_same_strand.txt"),
                                unsorted_reverse_counts = c("db/gene_assignment/ORFtag/Activator2_input_rep1_rev_strand.txt",
                                                            "db/gene_assignment/ORFtag/Activator2_input_rep2_rev_strand.txt"),
                                exons_start_gtf = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
                                name = "Activator_2",
                                output_suffix = "_vs_input_strandBias.txt",
                                output_folder_FC_file = "db/FC_tables/ORFtag/")
}
