#' Aligns sanger ab files to a refseq
#'
#' This function plots a heatmap showing alignment to a reference sequence
#'
#' @param refseq A character string corresponding to the reference sequence
#' @param abfiles Paths to abfiles to analyse
#' @param revcomp Boolean value specifying if the sequence shoul be reversed before alignment to the refseq
#' @param feat.mismatches Number of mismatches allowed to align features
#' @param feat.sequences additional feature sequences to plot on top of the matrix
#' @param feat.names Names of additional feature sequences
#' @param feat.cols Colors to use to plot additional feature sequences
#' 
#' @examples
#' refseq <- "CGCGCGCGTTGACAGTGAGCGCGTCTCTCACCGGAGCAAGTACTCCGTTCGAAGTTTTATACCGCTTAACTATTTAGAGGATGCCAGCCAACTGGATACGGATTTATTGACCAAATGCTGCAGTGACGTTCATCCGTGTGCCTTGCTCTTCTGCTACTTCAACAATCGAGGAGCTTTTCTGGACTATCTAAATGTGTATCAAGGAAAGTGGCTTATATTGATTGGACCATTGCCCGCCCTGGGAATCCACACAGATCCCAATCCCGCCCAGCCCCAACTGCCCCTAGGATCGACGCGGACAAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCCCGGGCCCGGGCCCCAGGCTTTCAACAAGTGCATTCTTAGTCCCCATATACCCTGCTTGTATGATTTTTTTGAACGAGGTATCGTTCAACTGGCAGAGAGTTCAAAATCATAAATTTCAAACGATATTAGGCCGCAGCCGCAGTGGCGCGTTGATGATGCAGCTCCGCTCCACTCTGCCCCACTGCCATTCTGCCACTCTGCCACCTGGTGGTGCCCCAAAGAGCAGACTGGTGCAGTTGCAGTTGCTGAGAATTTTTGGCCCATGCTGGTGCGTAGATTCCCCCCCTCACCCTTCCCATTGCACGGCTAGAGCGTGCAGCACATTTTCTTTTGTTGTACGAATGCCGGAGTGGATTTATTGATTTTGGTCAAAAGCCATGCTTAAGCGTAAATTATGTGCTGGCTTAAACGCTGCGCTCGAACTTTTCGCGTCTCTCTCGCTTGGCGGCAAGTGCGATGAAAGCCACTGATTTGTGTGCCAGTCTAGGGTACGGCCATAAATCCGTTCAATTTCATAACCAAACAGCAATTTATCAACCGCACCCCCCTTTGCGCAAAAACACCGCCCATATGCCAACGGACTAATTAAAAGCAGCGACAAATATGCGGCGCAAAGTGCGAAAATCGAAATGAAAATGAAAATGAAACGACGATGAAAAGCTAAGGGGCACAACAAAAAAAAGGCTGGGAAACTGTGTCAGAAAAGCCGGCCAAACGACAAGTAGCTGATATCGCATAAAGCAACGATTTAATTAGGCGCTGATAATTATTGTACAACTAGTTCTGGCCAACTTCGCCATCGTTCCGAGTTCAAAATCGAGTTTGGTCTCAGACTTTCAGGTTGAGGTAGGCGCAACCGACGATAGCCGCGGCCGCCTATCGCACATGTTAGATTTGGCTAATTTACGAAATTGGTAGTGCCCGATATCGCCGAGACCCGGCTAATTCAATCATTCGCTGCCGGCGTCTGGCATTAGCGGATGATTTATGGTAACATTTCGGACTCAAGCGAGCTCGATTCTGAGCTTAGTGCCTAGATTCTGGAATCAACATGCCAAGCTAACCACATAAATAGTATACAACCCATGCTATAAATAACACTGCATTTTAAAGATATTCAAAGTATAGAAGATAATGTTAATTATTTAAACCAACCTAAGTACAGATAAGTACTTAGGATAATTGTTATGTAAAAAAAAGTCATACCACCATTTATACAGAAATACACACCACTTATCGATTGAGTTACTATGTTATAAAACGAACTTTGTTTTAATACCGTGCCTATAACTCAAGCTTAATTTATATATTCCAGCTCGCATCGCCTGTCTCTACGGTTATTTTTTTTTTTTTTGAGCCGCCAGCGGCTTTCGTTGAAAAAACTTGGAACAAATTCGAATTGCTACCGCCAATCCATCGATTCTCGATCGAAGCTTAACCCAATTGATCTTACAACGCGAATGGCATGTTAAACACCCTACACTCGTGTAATTTGTGTCCGAAATGAAATAAATTTATGAACATCGATGGCAAAGTTGATTAGACGAAAGTCGCCGCCAACTCGGACATGGTCGTGGTCCATGGAGTCGCGGACATGGACCGGGATTTGGACTTGGATTTGGGCATGGCACGGTATGCTATGGTTGGGGGCCACCGAATAACCAGTTTGCCTCTGCTTAGCGGACACTTGTAATATTTCTGTTGAACTTCTCCTAATTGCAATTTGATGTTTCAGTTTTCGCTCTTTTGGTAATTAAGTTGTATTTGGATAAGGATTTTCTTTTTCCCTTTTTTATTTGAATCCATTTTTCGGGGCGCACTTAGCGGCAGCGTGTGAAAAGTTGTGACGGGGGCCCCCGGGGGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTTGACAGTGAGCGCGTCTCTCACCGATTCTCGCCGTTCAGCCGGCTTTTTTCTCAATTGTATCGTAATCGTGAGTCAGCATAAAGCCATCGCGCCTCTGTTGATCAACTCATTCTAGAAAGATACTCCAAACATATATGTACTATGCGATAGCTCTCCCCTTTCGATGATGATAAGAATGAGCATGAGTGGACTTCATTATGGGTCAGCACATTTGGTTTGGGCTGGGGGCTTCCAAAGAGGGGGAACACCACATTGTGACACACTCGTTTAAAGCTTTTGAAGCGTGCAGAAT"
#' 
#' abfiles <- c("/groups/stark/vloubiere/projects/pe_STARRSeq/db/sanger_sequencing/vllib_006_011_sanger/LOUB_vllib006.1_CASeq044_A01.ab1",
#' "/groups/stark/vloubiere/projects/pe_STARRSeq/db/sanger_sequencing/vllib_006_011_sanger/LOUB_vllib006.1_CASeq003_A04.ab1")
#' revcomp <- c(F, T)
#' feat.names <- c("illumina_F", "Flink_-2", "illumina_R", "R1link+0", "R2link+0", "R3link-3", "spacer")
#' feat.sequences <- c("ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
#' "GACAGTGAGCGCGTCTCTCACCG", 
#' "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", 
#' "TTGTCCGCGTCGATCCTAGG", 
#' "ATTCTGCACGCTTCAAAAGC", 
#' "AGCTGCTGTCTAAGGCACAGG", 
#' "CCGCGCGCGCTCATCAATGTATCTTATCATGTCTGCTCGAAGCGGCCGGCCGAATTCG")
#' feat.cols <- c("blue", "green", "red", "pink", "purple", "gold", "black")
#' vl_sanger_align(refseq = refseq,
#' abfiles = abfiles,
#' revcomp= c(F, T),
#' feat.names = feat.names,
#' feat.sequences = feat.sequences,
#' feat.cols= feat.cols)
#' 
#' @export
vl_sanger_align <- function(refseq, 
                            abfiles, 
                            revcomp= NULL,
                            feat.mismatches= 0, 
                            feat.sequences= NULL,
                            feat.names= NULL,
                            feat.cols= NULL)
{
  if(length(refseq)!=1)
    stop("length(refseq)!=1")
  if(any(is.na(abfiles)))
    stop("Some of the abfiles are NA!")
  if(any(!file.exists(abfiles)))
    stop("Some of the abfiles do not exist!")
  if(is.null(revcomp))
    revcomp <- rep(T, length(abfiles))
  if(length(revcomp) != length(abfiles))
    stop("revcomp vector length should match the length of provided abfiles")
  if(!is.null(feat.sequences))
  {
    if(any(is.na(feat.sequences)))
      stop("Some feat sequences are NAs -> stop")
    if(is.null(feat.names))
      feat.names <- paste0("feat_", seq(feat.sequences))
    if(length(feat.names) != length(feat.sequences))
      stop("length feat names should match the length of feat seqauences")
    if(is.null(feat.cols))
    {
      if(length(feat.sequences)==1)
        feat.cols <- "black" else
          feat.cols <- colorRamps::matlab.like2(length(feat.sequences))
    }
    if(length(feat.cols) != length(feat.sequences))
      stop("length feat colors should match the length of feat seqauences")
  }

  # Import sequences and check number of usable nt
  seq <- sapply(abfiles, function(x) sangerseqR::primarySeq(sangerseqR::readsangerseq(x)))
  keep <- sapply(seq, 
                 function(x) 
                   length(which(unlist(strsplit(as.character(x), "")) %in% c("A", "T", "C", "G")))>100)
  seq <- seq[keep]
  abfiles <- abfiles[keep]
  revcomp <- revcomp[keep]
  seq[revcomp] <- sapply(seq[revcomp], Biostrings::reverseComplement)
  seq <- Biostrings::DNAStringSet(c(refseq, sapply(seq, as.character)))
  
  # Perform multiple alignment
  require(msa)
  align <- msa::msa(seq, order = "input", method= "ClustalOmega") 
  align <- data.table::rbindlist(lapply(as.character(align@unmasked), function(x) strsplit(x, "")), idcol = T)
  align[, idx:= .SD[, .I], .id]
  align <- data.table::dcast(align, idx~.id, value.var = "V1")
  colnames(align)[-1] <- c("refseq", gsub(".ab1", "", basename(abfiles)))
  
  # make matrix
  mat <- as.matrix(align, 1)
  mat[mat=="-"] <- "white"
  for(i in 2:ncol(mat))
    mat[mat[, i] != "white", i] <- ifelse(mat[mat[, i] != "white", i]==mat[mat[, i] != "white", 1], "blue", "red")
  mat[mat[, 1]!="white",1] <- apply(mat[mat[, 1]!="white",], 1, function(x) ifelse(any(x=="blue"), "blue", "red"))
  
  # Plot figure ----
  plot.new()
  rasterImage(t(mat), 0, 0, 1, 1, interpolate = F)
  abline(h= seq(0, 1, length.out= ncol(mat)+1))
  box()
  y <- seq(1, 
           0, 
           length.out= ncol(mat)*2+1)
  y <- y[(seq(y) %% 2)==0]
  axis(side = 2, 
       at = y, 
       labels = colnames(mat), 
       las= 1)
  
  # Add features ----
  if(!is.null(feat.sequences))
  {
    feat <- data.table::data.table(name= feat.names,
                                   seq= feat.sequences,
                                   Cc= feat.cols)
    feat[, {
      # Add forward feature
      .c <- data.table::as.data.table(Biostrings::matchPattern(Biostrings::DNAString(seq),
                                                               subject = refseq, 
                                                               with.indels= T,
                                                               max.mismatch = feat.mismatches)@ranges)
      if(nrow(.c)>0)
      {
        arrows(.c$start/nrow(mat), 1.1, 
               .c$end/nrow(mat), 1.1, 
               lwd= 3, 
               col = Cc[1], 
               xpd= T, 
               lend= 2, 
               length = 0.025)
      }
      # Add reverse feature  
      .c <- data.table::as.data.table(Biostrings::matchPattern(Biostrings::reverseComplement(Biostrings::DNAString(seq)),
                                                               subject = refseq,
                                                               with.indels= T,
                                                               max.mismatch = feat.mismatches)@ranges)
      if(nrow(.c)>0)
      {
        arrows(.c$end/nrow(mat), 1.1, 
               .c$start/nrow(mat), 1.1, 
               lwd= 3, 
               col = Cc[1], 
               xpd= T, 
               lend= 2, 
               length = 0.025)
      }
    }, (feat)]
  }
  
  print("DONE")
}
