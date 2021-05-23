# require(data.table)
# require(ontologyIndex)
# 
# test <- readRDS("Rdata/final_clustering_transcriptomes.rds")
# test <- lapply(split(test, by= "cl"), function(x) x$FBgn)
# 
# vl_GO_cluster <- function(dat_list, 
#                           type,
#                           plot= T,
#                           padj_cutoff= 0.001,
#                           log2OR_cutoff_FUN= function(x) x>1,
#                           
#                           go_object= )
# {
#   if(!is.list(dat_list))
#     stop("dat_list must be a list object") else
#       dat <- rbindlist(lapply(dat_list, as.data.table), idcol = "cl")
#   if(!type %in% c("FBgn", "Symbol"))
#     stop("Type must be either FBgn or Symbol")
#   dat[, cl:= as.factor(cl)]
#   colnames(dat)[2] <- type
#   
#   # Compute enrichment
#   go <- copy(go_object)
#   go[, genes_GO:= length(unique(Fbgn)), GO]
#   go[, genes_noGO:= length(unique(go$Fbgn))-genes_GO]
#   enrich <- dat[, {
#     .c <- unique(FBgn)
#     res <- go[.c, .(genes_GO, 
#                     genes_noGO, 
#                     cl_GO= length(unique(Fbgn))), 
#               on= "Fbgn", 
#               by= .(GO, name), 
#               nomatch= NULL]
#     res[, cl_noGO:= length(.c)-cl_GO]
#     res <- unique(res)
#     res[, fisher.test(data.table(c(cl_GO, cl_noGO), 
#                                  c(genes_GO, genes_noGO)))[c("estimate", "p.value")], (res)]
#   }, cl]
#   enrich[, p.adj:= p.adjust(p.value, method = "fdr"), cl]
#   enrich[, log2OR:= log2(estimate)]
#   
#   # Select signif ontologies
#   sel <- enrich[, any(p.adj<=padj_cutoff & log2OR_cutoff_FUN(log2OR)), GO][(V1), GO]
#   dmat <- enrich[GO %in% sel]
#   setorderv(dmat, c("cl", "log2OR"), c(1, -1))
#   dmat[, x:= .GRP, .(cl, log2OR)]
#   dmat[, x:= min(.SD[p.adj<=padj_cutoff, x]), name]
#   dmat[, x:= .GRP, keyby= x]
#   setorderv(dmat, "cl", order = -1)
#   dmat[, y:= .GRP, cl]
#   
#   # PLOT
#   Cc <- cut(dmat$p.adj, 
#             c(-Inf, 0.00001, 0.001, 0.01, 0.05, Inf), 
#             c("red", "tomato", "cornflowerblue", "blue", "black"))
#   
#   
#   plot.new()
#   .v <- seq(0, 1, length.out = max(dmat$x))
#   segments(.v, 0, .v, 1, lwd= 0.5)
#   .h <- seq(0, 1, length.out = max(dmat$y))
#   segments(0, .h, 1, .h, lwd= 0.5)
# 
#   points(x= .v[dmat$x], 
#          y= .h[dmat$y], 
#          cex= dmat$log2OR/2.5, 
#          col= as.character(Cc), 
#          pch= 16, 
#          xpd= T)
#   axis(2, 
#        at = .h, 
#        labels = unique(dmat$cl), 
#        las= 1, 
#        tick= 0, 
#        lwd= 0)
#   par(cex= 0.7)
#   axis(1, 
#        at = .v, 
#        labels = unique(dmat[, .(name, x)])$name, 
#        las= 2, 
#        tick= 0, 
#        lwd= 0)
# }