#required objects

#results_dge
#dge_rseq_fit
#p6c_data
#p6d_data
#p6e_data
#sAnalysisName 
#sSign
#lScoreENTREZ
#results_tfs


library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(ggplotify)

load("figure5_objects.rda")


######################################################
# Figure 5A
######################################################
p5a_data <- results_dge[["logcpm-wb-0.05_limma_ABZB_KO - ABZB_SCR"]]
p5a_data[,"label"] <- factor(ifelse(
  p5a_data[,"selected"] & p5a_data[,"logFC"] < 0,
  "DOWN",
  ifelse(
    p5a_data[,"selected"] & p5a_data[,"logFC"] > 0,
    "UP",
    "Not d.e."
  )),
  levels = c("UP", "DOWN", "Not d.e.")
)
rownames(p5a_data) <- p5a_data[,"Row.names"]
dataset_label <- "logcpm-wb-0.05_limma"
p5a <- as.ggplot(as.grob(function() {
  limma::plotMA(
    dge_rseq_fit[[dataset_label]],
    coef = 1, 
    status= p5a_data[sapply(
      rownames(dge_rseq_fit[[dataset_label]]$coefficients),
      function(s) strsplit(s,".",fixed=TRUE)[[1]][1]),
      "label"
    ],
    hl.col=c("red","blue","gray")
  )
  selected_genes <- c(
    head(results_dge[[paste0(dataset_label,"_ABZB_KO - ABZB_SCR")]][,"Row.names"],10),
    "ENSG00000178607"
  )
  names(selected_genes) <- results_dge[[paste0(dataset_label,"_ABZB_KO - ABZB_SCR")]][
    results_dge[[paste0(dataset_label,"_ABZB_KO - ABZB_SCR")]][,"Row.names"] %in% selected_genes,
    "gene_symbol"
  ]
  names(selected_genes)[is.na(names(selected_genes))] <- selected_genes[is.na(names(selected_genes))]
  selected_genes <- sapply(
    selected_genes,
    function(s) grep(
      s,
      names(dge_rseq_fit[[dataset_label]]$Amean),
      value = TRUE,
      fixed=TRUE
    )
  )
  text(
    x = dge_rseq_fit[[dataset_label]]$Amean[selected_genes],
    y = dge_rseq_fit[[dataset_label]]$coefficients[selected_genes,"ABZB_KO - ABZB_SCR"],
    labels = names(selected_genes)
  )
}))

pdf("Figure_5A.pdf")
print(p5a)
dev.off()





######################################################
# Figure 5B
######################################################

# dataset_label <- "cpmrbe_AB_all"
# dataset_label <- "vst_AB_all"
# dataset_label <- "logcpm_AB_all"
# p5b_data <- cbind(
#   cbind(
#     results_tfs[[dataset_label]][,"TF",drop=FALSE],
#     results_tfs[[dataset_label]][,c("S1B_SCR","S2B_SCR","S1B_G8_5KO","S2B_G8_5KO")]
#     # apply(
#     #   results_tfs[[dataset_label]][,c("S1B_SCR","S2B_SCR","S1B_G8_5KO","S2B_G8_5KO")],
#     #   2,
#     #   function(v) (v-min(v))/(max(v)-min(v))
#     # )
#   ),
#   results_tfs[[dataset_label]][,"ABZB_KO - ABZB_SCR p.val",drop=FALSE]
# )
# p5b_data[,"ABZB_KO - ABZB_SCR p.val"] <- unlist(apply(
#   p5b_data[,c('S1B_G8_5KO','S2B_G8_5KO','S1B_SCR','S2B_SCR')],
#   1,
#   function(o) {
#     ret.p.val <- t.test(
#       o[c('S1B_G8_5KO','S2B_G8_5KO')],
#       o[c('S1B_SCR','S2B_SCR')],
#       paired = TRUE
#     )$p.value
#     return(ret.p.val)
#   } 
# ))

# p5b <- ggplot(p5b_data) + geom_point(
#   aes(
#     x = ((S1B_G8_5KO + S2B_G8_5KO)) / (S1B_SCR + S2B_SCR),
#     y = -log10(get("ABZB_KO - ABZB_SCR p.val")),
#     fill= ifelse(get("ABZB_KO - ABZB_SCR p.val")<0.05,get("ABZB_KO - ABZB_SCR p.val"),NA)
#   ),
#   shape = 21,
#   size = 4
# ) + geom_label_repel(
#   aes(
#     x = ((S1B_SCR + S2B_SCR) - (S1B_G8_5KO + S2B_G8_5KO)) / (S1B_SCR + S2B_SCR),
#     y = -log10(get("ABZB_KO - ABZB_SCR p.val")),
#     label= ifelse(get("ABZB_KO - ABZB_SCR p.val") <= 0.05,TF,NA)
#   )
# ) + geom_label_repel(
#   aes(
#     x = ((S1B_SCR + S2B_SCR) - (S1B_G8_5KO + S2B_G8_5KO)) / (S1B_SCR + S2B_SCR),
#     y = -log10(get("ABZB_KO - ABZB_SCR p.val")),
#     label= ifelse(TF == "ATF6", TF, NA)
#   ),
#   color = "red"
# ) + scale_fill_gradient(
#   name = "p-val.",
#   trans = "log",
#   low = "red",
#   high = "white",
#   labels = c("0.05","0.01","0.001"),
#   breaks = c(0.05,0.01,0.001)
# ) + scale_color_discrete(guide="none") + xlab(
#   "Fraction of change"
# ) + ylab(
#   "p-value"
# ) + theme_bw() + theme(
#   axis.text.x = element_text(size = 16),
#   axis.text.y = element_text(size = 16),
#   axis.title.x =element_text(size = 16),
#   axis.title.y = element_text(size = 16),
#   legend.key.size = unit(0.5, 'cm'),
#   legend.text = element_text(size=16)
# )
# plot(p5b)



######################################################
# Figure 5C
######################################################
sCategory <- "logcpm-wb-0.05_limma_ABZB_KO - ABZB_SCR_UP_enrichGOBP"
p5c <- dotplot(p6c_data[,asis=T], orderBy = "x", showCategory=10) + ggtitle(
      paste("Dotplot for",sAnalysisName,sSign,sCategory,'for enrich sign',"UP"))

pdf("Figure_5C.pdf")
print(p5c)
dev.off()



######################################################
# Figure 5D
######################################################
sCategory <- "logcpm-wb-0.05_limma_ABZB_KO - ABZB_SCR_DOWN_enrichGOBP"
    
p5d <- dotplot(p6d_data[,asis=T], orderBy = "x", showCategory=10) + ggtitle(
      paste("Dotplot for",sAnalysisName,sSign,sCategory,'for enrich sign',"DOWN"))
 
pdf("Figure_5D.pdf")
print(p5d)
dev.off() 


######################################################
# Figure 5E
######################################################



p5e <- as.ggplot(cnetplot(
      p6e_data, 
      showCategory=10, 
      categorySize="pvalue", 
      foldChange = lScoreENTREZ, 
      colorEdge = TRUE
    ))

pdf("Figure_5E.pdf")
print(p5e)
dev.off() 



