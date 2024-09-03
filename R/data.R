#' Dmel motifs data base 
#'
#' Contains 13899 motifs usable for motif enrichments in Dmel
#'
#' @usage Can be used with the ?vl_motif_counts() function to count motifs
#' 
#' @format An object containing 13899 motifs and related metadata
#' \describe{
#'   \item{motif_name}{Motif ID from TF_clusters_PWMs.RData}
#'   \item{FBgn}{Curated FBgn symbols. sep= "/"}
#'   \item{motif_cluster}{Motif cluster from Bernardo Almeida}
#'   \item{collection}{}
#'   \item{collection_version}{}
#'   \item{species}{PWM expressed as log odds}
#'   \item{pwms_log_odds}{list of PWM matrices}
#'   \item{pwms_perc}{list of PWM matrices}
#' }
#' @source {"/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData"}
"vl_Dmel_motifs_DB_full"

#' Motifs DB noin-redundant version (BA) 
#'
#' Contains 6502 non-redundant motifs
#'
#' @usage This curated list of non-redundant motifs was assembled by Bernardo P. de almeida (https://github.com/bernardo-de-almeida/motif-clustering). It can be used with the ?vl_motif_counts() function to count motifs
#' 
#' @format An object containing 6502 motifs and related metadata
#' \describe{
#'   \item{motif_name}{Motif ID from TF_clusters_PWMs.RData}
#'   \item{motif_cluster}{Motif cluster from Bernardo Almeida}
#'   \item{collection}{}
#'   \item{collection_version}{}
#'   \item{species}{PWM expressed as log odds}
#'   \item{pwms_log_odds}{list of PWM matrices}
#'   \item{pwms_perc}{list of PWM matrices}
#' }
#' @source {"/groups/stark/almeida/Papers/DeepSTARR/Code/TF_motif_database/TF_clusters_PWMs.RData"}
"vl_motifs_DB_v2"

#' SUHW top peaks
#'
#' Example set that can be used for motif enrichment...
#'
#' @usage see ?vl_motif_enrich()
#' @format Narrowpeaks file containing top SUHW peaks on dm3 chanonical chromosomes
"vl_SUHW_top_peaks"

#' STARR-Seq top peaks
#'
#' Example set that can be used for motif enrichment...
#'
#' @usage see ?vl_motif_enrich()
#' @format Narrowpeaks file containing top STARR-Seq peaks (DSCP) on dm3
"vl_STARR_DSCP_top_peaks"

#' Set of genes
#'
#' Example set that can be used for GO enrichment. Contain RpL and HOX genes
#'
#' @usage see ?vl_GO_enrich()
#' @format FBgn, symbol and GO (RpL/HOX)
"vl_genes_set"

#' Thermofisher enzymes
#'
#' Contains cutting sites for most commercial enzymes
#'
#' @usage Can be used with vl_digest. see ?vl_digest()
#' 
#' @format A table containing thermofisher enzymes and related info
#' \describe{
#'   \item{name1}{enzyme name}
#'   \item{name2}{Alternative names}
#'   \item{fastdigest}{FastDigest?}
#'   \item{buffer}{Buffer to be used}
#'   \item{temperature}{Temeprature for digest}
#'   \item{cutsite}{cutsite pattern}
#'   \item{consensus_F}{consensus motif}
#'   ...
#' }
#' @source {"https://www.thermofisher.com/at/en/home/brands/thermo-scientific/molecular-biology/thermo-scientific-restriction-modifying-enzymes/restriction-enzymes-thermo-scientific/conventional-restriction-enzymes-thermo-scientific.html"}
"vl_thermofisher_restriction_enzymes_table"

#' Example metadata ORFtag
#'
#' Metadata table example for the ORFtag pipeline. See ?ORFtag_pipeline() for further details
#'
#' @usage See ?ORFtag_pipeline()
#' 
#' @format data.table
#' #' \describe{
#'   \item{user}{User name}
#'   \item{batch}{batch}
#'   \item{screen}{screen name}
#'   \item{condition}{Typically one of 'sort' or 'input'}
#'   \item{replicate}{replicate}
#'   \item{barcodes}{Typically two barcodes, e.g. 'GCCTCTTC|ATTGATTC'}
#'   \item{sampleID}{sampleID. Should correspond to the catenation of screen, condtion and replicate  (separated by '_'), and be unique for each biological sample/replicate. Samples with the same sampleID will be merged (re-sequencing...).}
#'   \item{expName}{Experiment name}
#'   \item{genome}{For now, only 'mm10' and 'hg38' are supported}
#'   \item{layout}{'PAIRED' or 'SINGLE'. Will be used in combination with 'sequencer' to know in which column barcodes should be found}
#'   \item{sequencer}{'NextSeq or HiSeq. Will be used in combination with 'layout' to know in which column barcodes should be found}
#'   \item{bam_path}{The path to the VBC bam file containing the reads}
#'   \item{comment}{Extra columns are tolerated}
#'   ...
#' }
"vl_metadata_ORFtag"

#' Example metadata CUTNRUN
#'
#' Metadata table example for the CUTNRUN pipeline. See ?vl_CUTNRUN_processing() for further details
#'
#' @usage See ?vl_CUTNRUN_processing()
#' 
#' @format data.table
#' #' \describe{
#'   \item{user}{User name}
#'   \item{method}{CutNRun, ChIP, ....}
#'   \item{experiment}{exp1, exp2, ...}
#'   \item{antibody}{V5, H3K27Ac, ...}
#'   \item{cell_line}{AID.Hcfc1, parental, ...}
#'   \item{treatment}{noTreatment, IAA, ...}
#'   \item{replicate}{rep1, rep2, ...}
#'   \item{condition}{Should be the concatenation of ab, cell_line, treatment and experiment  (separated by '_'). Therefore, is is similar to sampleID but without the replicate information, and will be used to call peaks on the merge data.}
#'   \item{sampleID}{sampleID should be the concatenation of ab, cell_line, treatment, experiment and replicate  (separated by '_'), and be unique for each sample/replicate. Samples with the same sampleID will be merged (re-sequencing...).}
#'   \item{barcode}{experimental barcode. e.g. ACAGTG, GCCAAT...}
#'   \item{i5}{i5 index, such as ACGTCCTG...}
#'   \item{genome}{For now, only 'mm10' and 'hg38' are supported}
#'   \item{input}{Sample ID of the input data to be used for callling peaks}
#'   \item{layout}{'PAIRED' or 'SINGLE'. Will be used in combination with 'sequencer' to know in which column barcodes should be found}
#'   \item{sequencer}{'NextSeq or HiSeq. Will be used in combination with 'layout' to know in which column barcodes should be found}
#'   \item{bam_path}{The path to the VBC bam file containing the reads}
#'   \item{comment}{Extra columns are tolerated}
#'   ...
#' }
"vl_metadata_CUTNRUN"

#' Example metadata PROseq
#'
#' Metadata table example for the CUTNRUN pipeline. See ?vl_PROseq_processing() for further details
#'
#' @usage See ?vl_CUTNRUN_processing()
#' 
#' @format data.table
#' #' \describe{
#'   \item{user}{User name}
#'   \item{experiment}{HCFC1, HCFC1_rescue, ...}
#'   \item{cell_line}{AID-Hcfc1-cl4, AID-Hcfc1-cl17, ...}
#'   \item{treatment}{0hrIAA, 3hrIAA, ...}
#'   \item{replicate}{rep1, rep2, ...}
#'   \item{sampleID}{sampleID should be the concatenation of cell_line, treament and replicate (separated by '_'), and be unique for each sample/replicate. Samples with the same sampleID will be merged (re-sequencing...).}
#'   \item{DESeq2_name}{Simplified sample name to be used within the DESeq2 object. Should be unique for each sample/replicate.}
#'   \item{DESeq2_condition}{Condition name to be used in DESeq2. Ideally, it will be similar to DESeq2_name but with no replicate information.}
#'   \item{DESeq2_control}{Should correspond to the DESeq2_condition of the control condition to be used for differential analysis.}
#'   \item{barcode}{Index barcode. e.g. ACAGTG, ...}
#'   \item{eBC}{Sample barcode found at the beginning of the read. e.g. ATCG, GCTA...}
#'   \item{genome}{For now, only 'mm10' and 'hg38' are supported}
#'   \item{spikein_genome}{For now, only 'dm3' is supported}
#'   \item{layout}{'PAIRED' or 'SINGLE'. Will be used in combination with 'sequencer' to know in which column barcodes should be found}
#'   \item{sequencer}{'NextSeq or HiSeq. Will be used in combination with 'layout' to know in which column barcodes should be found}
#'   \item{bam_path}{The path to the VBC bam file containing the reads}
#'   \item{comment}{Extra columns are tolerated}
#'   ...
#' }
"vl_metadata_PROseq"