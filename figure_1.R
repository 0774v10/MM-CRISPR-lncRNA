#require objects: 
#CR_disp
#table_MAGeCK_sgrna
#gstable_ABZB
#gstable_AMO
#gene_intersections

library(ggplotify)
library(pheatmap)
library(ggplot2)
load("figure1_objects.rda")

############################################################
# figure 1B
############################################################

p1b <- as.ggplot(pheatmap(
      CR_disp[c(-1,-2),c(-1,-2)],
      cluster_cols = TRUE,
      show_colnames = TRUE,
      filename = 'Figure_1B.pdf',
      color = colorRampPalette(c('black','green'))(100),
      annotation_col = data.frame("CL"=sapply(colnames(CR_disp)[-2:-1],function(s) strsplit(s,"_")[[1]][1])),
      annotation_row = data.frame("CL"=sapply(rownames(CR_disp)[-2:-1],function(s) strsplit(s,"_")[[1]][1])),
      annotation_names_col = FALSE,
      annotation_names_row = FALSE,
      # labels_col = sapply(colnames(CR_disp)[-2:-1],function(s) strsplit(s,"_")[[1]][2]),
      # labels_row = sapply(rownames(CR_disp)[-2:-1],function(s) strsplit(s,"_")[[1]][2]),
      labels_row = paste0(c(rep("AMO1 r",3),rep("ABZB r",3)),c(1:3)),
      labels_col= paste0(c(rep("AMO1 r",3),rep("ABZB r",3)),c(1:3)),
      annotation_colors = list(CL=c(AMO="blue",ABZB="goldenrod1")),
      angle_col = 45,
      # legend = FALSE,
      annotation_legend = FALSE,
      display_numbers = FALSE,
      cellwidth = 40,
      cellheight = 40,
      fontsize = 20,
      treeheight_row = 0, 
      treeheight_col = 0
))

pdf("Figure_1B.pdf")
plot(p1b)
dev.off()

############################################################
# figure 1C
############################################################



p1c <- ggplot(table_MAGeCK_sgrna[
  !is.na(table_MAGeCK_sgrna$Type) &
  !is.na(table_MAGeCK_sgrna$SubType) &
  table_MAGeCK_sgrna$Type!="AAVS1" &
  table_MAGeCK_sgrna$SubType!="ncRNAs introns" &
  table_MAGeCK_sgrna$SubType!="Ribosomal Proteins\nintrons",
])  + geom_jitter(
aes(y=LFC, x=SubType), #fill=SubType),
  color = "blue4",
  alpha=0.5,
  shape = 19,
  size = 2
) + facet_grid(
  cell_line ~ .,
  scales = "free_x"
) + geom_boxplot(color="grey60",fill="white",lwd=1,
  aes(y=LFC, x=SubType, color=SubType),
  outlier.shape = NA
)+ theme_bw() + guides(
  color="none",
  fill="none"
) + labs(x = NULL) + theme(
  strip.text = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 20),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 20)
)
pdf("Figure_1C.pdf")
plot(p1c)
dev.off()


############################################################
# figure 1D
############################################################

nThresholdFDR=0.2

#rank gstable_ABZB
gstable_ABZB_rnk=gstable_ABZB[order(gstable_ABZB$ensembl_transcript_id),]
#rank gstable_AMO
gstable_AMO_rnk=gstable_AMO[order(gstable_AMO$ensembl_transcript_id),]
#rank gene_intersections
gene_intersections_rnk=gene_intersections[order(gene_intersections$ensembl_transcript_id.x),]

#join the missing column (gene_biotype) from gene_intersections_rnk to gstable_ABZB_rnk and gstable_AMO_rnk
gstable_AMO_rnk$gene_biotype=gene_intersections_rnk$gene_biotype
gstable_ABZB_rnk$gene_biotype=gene_intersections_rnk$gene_biotype

#use gstable_AMO_rnk and gstable_ABZB_rnk to produce the figure 1D
#select not controls, not ribosomal and with neg.fdr<=nThresholdFDR for each ABZB and AMO1
#also rerank according to neg rank
ABZB_toplot_logical= !gstable_ABZB_rnk$is_control &!gstable_ABZB_rnk$is_ribosomal&!gstable_ABZB_rnk$gene_biotype=="protein_coding" & gstable_ABZB_rnk$neg.fdr  <=nThresholdFDR
AMO1_toplot_logical= !gstable_AMO_rnk$is_control &!gstable_AMO_rnk$is_ribosomal&!gstable_AMO_rnk$gene_biotype=="protein_coding" & gstable_AMO_rnk$neg.fdr  <=nThresholdFDR

significance=ifelse(ABZB_toplot_logical &AMO1_toplot_logical,"both",ifelse(ABZB_toplot_logical & !AMO1_toplot_logical, "ABZBonly",ifelse(!ABZB_toplot_logical & AMO1_toplot_logical,"AMO1only","none")))
gstable_ABZB_rnk$hit_type=significance
gstable_AMO_rnk$hit_type=significance

gstable_ABZB_rnk=gstable_ABZB_rnk[order(-gstable_ABZB_rnk$neg.score),]
gstable_AMO_rnk=gstable_AMO_rnk[order(-gstable_AMO_rnk$neg.score),]

pdf("Figure_1D.pdf",width=5,height=10)
par(mfrow=c(2,1))
plot(x=gstable_ABZB_rnk$neg.rank,y=-log10(gstable_ABZB_rnk$neg.score),ylab="-log10 MAGeCK RRA score",xlab="MAGeCK RRA rank",
      type="l",lwd=2,col="grey50",main="ABZB hits",xlim=c(0,60))
points (x=gstable_ABZB_rnk$neg.rank[gstable_ABZB_rnk$hit_type=="ABZBonly"],y=-log10(gstable_ABZB_rnk$neg.score)[gstable_ABZB_rnk$hit_type=="ABZBonly"],type="p",pch=19,
      col="red")
points (x=gstable_ABZB_rnk$neg.rank[gstable_ABZB_rnk$hit_type=="both"],y=-log10(gstable_ABZB_rnk$neg.score)[gstable_ABZB_rnk$hit_type=="both"],type="p",pch=19,
      col="purple")
#ordercol=order(p1d_data_toplot[!is.na(p1d_data_toplot$Top),]$neg.lfc)
legend("topright",pch=19,legend=c("ABZB only","ABZB and AMO-1"),col=c("red","purple"))

plot(x=gstable_AMO_rnk$neg.rank,y=-log10(gstable_AMO_rnk$neg.score),ylab="-log10 MAGeCK RRA score",xlab="MAGeCK RRA rank",
      type="l",lwd=2,col="grey50",main="AMO-1 hits",xlim=c(0,30))
points (x=gstable_AMO_rnk$neg.rank[gstable_AMO_rnk$hit_type=="AMO1only"],y=-log10(gstable_AMO_rnk$neg.score)[gstable_AMO_rnk$hit_type=="AMO1only"],type="p",pch=19,
      col="steelblue1")
points (x=gstable_AMO_rnk$neg.rank[gstable_AMO_rnk$hit_type=="both"],y=-log10(gstable_AMO_rnk$neg.score)[gstable_AMO_rnk$hit_type=="both"],type="p",pch=19,
      col="purple")
#ordercol=order(p1d_data_toplot[!is.na(p1d_data_toplot$Top),]$neg.lfc)
legend("topright",pch=19,legend=c("AMO-1 only","AMO-1 and ABZB"),col=c("steelblue1","purple"))
dev.off()
################################################################################















############################################################
# figure 1E and Suppl. table S6
############################################################


inputList <- list(
ABZB = gstable_ABZB_rnk[
  !gstable_ABZB_rnk[,'is_control'] & !gstable_ABZB_rnk$gene_biotype=="protein_coding" &
     !gstable_ABZB_rnk[,'is_ribosomal'] &
    gstable_ABZB_rnk[,'neg.fdr'] <= nThresholdFDR,
  'id'],
"AMO-1" = gstable_AMO_rnk[
  !gstable_AMO_rnk[,'is_control'] & !gstable_AMO_rnk$gene_biotype=="protein_coding" &
     !gstable_AMO_rnk[,'is_ribosomal'] & 
    gstable_AMO_rnk[,'neg.fdr'] <= nThresholdFDR & gstable_AMO_rnk$id!="AAVS1" & supp_tab_AMO1$id!="HOTAIR",
  'id']
)

p1e <- ggvenn::ggvenn(
      inputList,stroke_size = 1,
      set_name_size = 12,
      show_percentage = FALSE,
      text_size = 12,
      fill_color=c("red","steelblue1")
) + theme(
      plot.margin = margin(t = 4, r = 0, b = 0, l = 0, unit = "pt")
)
pdf("Figure_1E.pdf")
plot(p1e)
dev.off()


