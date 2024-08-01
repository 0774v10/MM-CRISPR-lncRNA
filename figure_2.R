#require objects: 
#p2c_data
#p2d_plotdata_details
#df_fig2F
#gene_intersections
#data_visits
#response_status
#data_survival
#gene_list_survival
#table_MAGeCK_sgrna

library(survival)
library(survminer)
library(dplyr)
library(ggnewscale)
library(ggrepel)

load("figure2_objects.rda")





############################################################
# figure 2B
############################################################
gene_id <- "ENSG00000228013"

p2b_data <- table_MAGeCK_sgrna[
  table_MAGeCK_sgrna[,"Gene"] %in% c(
    gene_intersections[
      !is.na(gene_intersections[,"gene_id"]) &
        gene_intersections[,"gene_id"] == gene_id,
      "id"
    ],
    "NonTargetControl"
  ),
]
p2b_data[,"Gene"] <- ifelse(
  p2b_data[,"Gene"] == "NonTargetControl",
  "Non-Target",
  sapply(p2b_data[,"Gene"],function(s) strsplit(s,"_",fixed=TRUE)[[1]][2])
)
p2b <- ggplot(p2b_data) + geom_boxplot(
  aes(y=LFC, x=SubType),
  color = "black",
  outlier.shape = NA
) + geom_jitter(
  aes(y=LFC, x=SubType, fill=SelectedProbe),
  color = "black",
  alpha=0.5,
  shape = 21,
  size = 4
) + facet_grid(
  cell_line ~ Gene, scales = "free_x"
) + scale_color_discrete(
  guide="none"
)  + theme_bw() + theme(
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  strip.text = element_text(size = 16),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)+scale_fill_manual(values = ifelse(p2b_data$SelectedProbe , "blue3", "red"))


pdf("Figure_2B.pdf")
plot(p2b)
dev.off()











############################################################
# figure 2C
############################################################


response_status <- "Baseline"    
# p2c_data <- t(data_expression_patients[
#   data_expression_patients[,'GENE_ID'] %in% unique(gene_intersections[,"gene_id"]) &
#   data_expression_patients[,'GENE_ID'] %in% gene_list_survival,
#   data_visits[data_visits[,'response_status']==response_status,'sample_id']
# ])
p2c_data <- cbind(
  data_survival[
    sapply(rownames(p2c_data),function(s) paste(strsplit(s,'_')[[1]][1:2],sep='_',collapse='_')),
    c('public_id','ttcpfs','censpfs','ttcos','censos')],
  p2c_data
)

# Run test for all the screened genes
p2c_coxph_all <- NULL
for (surv_type in c("os","pfs")) {    
  for (gene_id in na.omit(unique(gene_intersections[,"gene_id"]))) {
    if (gene_id %in% colnames(p2c_data)) {
      p2c_coxph <- suppressWarnings(coxph(formula(paste0('Surv(ttc',surv_type,', cens',surv_type,') ~ ',gene_id)), data=p2c_data))
      p2c_coxph_all <- rbind(p2c_coxph_all, cbind(
        data.frame(
          gene_id = gene_id,
          surv_type = surv_type,
          coef = summary(p2c_coxph)[["coefficients"]][,"coef"],
          exp_coef = summary(p2c_coxph)[["coefficients"]][,"exp(coef)"],
          test_type  = c('waldtest','logtest','sctest'),
          stringsAsFactors = FALSE
        ),
        do.call(rbind,summary(p2c_coxph)[c('waldtest','logtest','sctest')])
      ))
    }
  }
}


# Set parameters of interest

# Estimate surv for the the gene of interest
print(paste("Survival", surv_type))
print(do.call(rbind,summary(p2c_coxph)[c('waldtest','logtest','sctest')]))
p2c_data <- p2c_data[order(p2c_data[,gene_id], decreasing = FALSE),]

# Get best split point for Kaplan-Maier curves
for (surv_type in c("os","pfs")) {
  dPval <- data.frame()
  for (sThr in c("median","mean",1:99)) {
    if (sThr %in% c("mean","median")) {
      if (sThr=="mean") nThr <- mean(p2c_data[,gene_id])
      if (sThr=="median") nThr <- median(p2c_data[,gene_id])
    } else {
      nThr <- quantile(p2c_data[,gene_id],as.integer(sThr) / 100)
    }
    nPval <- surv_pvalue(survfit(formula(paste0('Surv(ttc',surv_type,', cens',surv_type,') ~ ',gene_id,' <',nThr)) , data=p2c_data))$pval
    dPval <- bind_rows(
      dPval,
      data.frame(gene_id = gene_id, threshold_type = sThr, threshold_val = nThr, pval = nPval)
    )
  }
  dPval[!dPval[,"threshold_type"] %in%c("mean","median"),"pval_ma"] <- forecast::ma(
    x = dPval[!dPval[,"threshold_type"] %in%c("mean","median"),"pval"],
    order = 3,
    centre = TRUE
  )
  dPval[dPval[,"threshold_type"] == "mean","pval_ma"] <- dPval[dPval[,"threshold_type"] == "mean","pval"]
  dPval[dPval[,"threshold_type"] == "median","pval_ma"] <- dPval[dPval[,"threshold_type"] == "median","pval"]
  dPval[,"selected"] <- !is.na(dPval[,"pval"]) & dPval[,"pval"] == min(dPval[,"pval"],na.rm= TRUE)
  print(paste("Best split point for",surv_type))
  print(dPval[dPval[,"selected"],])
  # write.table(
  #   dPval,
  #   file = paste0("./figures/Figure_2a_surfit_",surv_type,".txt"),
  #   row.names = FALSE,
  #   col.names = TRUE,
  #   sep = "\t"
  # )
}


p2c <- ggarrange(
  ggsurvplot(
    survfit(formula(paste0('Surv(ttc',"os",', cens',"os",') ~ ',gene_id,' <',0.2784258)) , data=p2c_data),
    pval = TRUE, 
    conf.int = TRUE,
    risk.table = TRUE, 
    risk.table.col = "absolute",
    legend.labs = c("Higer","Lower"),
    linetype = "strata", 
    surv.median.line = "hv",
    ggtheme = theme_bw(), 
    palette = c("#E7B800", "#2E9FDF")
  )$plot + ylab("Survival prob. (OS)") + theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x =element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.key.size = unit(0.1, 'cm'),
    legend.text = element_text(size=16)
  ),
  ggsurvplot(
    survfit(formula(paste0('Surv(ttc',"pfs",', cens',"pfs",') ~ ',gene_id,' <',0.227183)) , data=p2c_data),
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "absolute",
    legend.labs = c("Higer","Lower"),
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_bw(),
    palette = c("#E7B800", "#2E9FDF")
  )$plot + ylab("Survival prob. (PFS)") + theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x =element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.key.size = unit(0.1, 'cm'),
    legend.text = element_text(size=16)
  ),
  # ggsurvplot(
  #   survfit(formula(paste0('Surv(ttc',"os",', cens',"os",') ~ ',gene_id,' <',0.2784258)) , data=p2c_data),
  #   pval = TRUE, 
  #   conf.int = TRUE,
  #   risk.table = 'absolute', 
  #   risk.table.col = "strata",
  #   legend.labs = c("Higer","Lower"),
  #   linetype = "strata", 
  #   surv.median.line = "hv",
  #   ggtheme = theme_bw(), 
  #   palette = c("#E7B800", "#2E9FDF")
  # )$table + scale_y_discrete(
  #   labels=c("Low ", "High")
  # ) + theme(
  #   axis.text = element_text(size = 16),
  #   axis.title = element_text(size = 16),
  #   legend.position="none"
  # ),
  ncol = 1
)



pdf("Figure_2C.pdf")
plot(p2c)
dev.off()




###########################################################
# FIGURE 2D
############################################################
# gene_id <- "ENSG00000228013" #(RP11-350G8.5)
# p2d_plotdata_details <- do.call(rbind,lapply(
# c("bone_marrow", "normal_tissue","patient_Baseline_low", "patient_Baseline_high"),
# function(sCol) cbind(
#   data.frame(
#     sample_type = switch(
#       sCol,
#       bone_marrow = "Bone Marrow",
#       normal_tissue = "Normal Tissue",
#       patient_Baseline_low = "Baseline (low)",
#       patient_Baseline_high = "Baseline (high)"
#     ),
#     stringsAsFactors = FALSE
#   ),
#   t(de_data[
#     gene_id,
#     rownames(design.matrix)[design.matrix[, gsub("_high", "", gsub("_low", "", sCol, fixed = TRUE), fixed = TRUE)]==1],
#     drop = FALSE
#   ])[
#     ifelse(
#       t(de_data[
#         gene_id,
#         rownames(design.matrix)[design.matrix[, gsub("_high", "", gsub("_low", "", sCol, fixed = TRUE), fixed = TRUE)]==1],
#         drop = FALSE
#       ]) >= 0.2784258,
#       if (sCol %in% c("bone_marrow", "normal_tissue")) {
#         TRUE
#       } else {
#         if (sCol == "patient_Baseline_high") {TRUE} else {FALSE}
#       },
#       if (sCol %in% c("bone_marrow", "normal_tissue")) {
#         TRUE
#       } else {
#         if (sCol == "patient_Baseline_high") {FALSE} else {TRUE}
#       }
#     ),, drop = FALSE
#   ]
# )))




#add healthy B cells 
gene_id <- "ENSG00000228013" #(RP11-350G8.5)
#extract only the 3 replicates of plasma cells (gene_id <- "ENSG00000228013", IL6R-AS1, alias RP11-350G8.5)
bulk_RNAseq_GSE148924=read.table("RNAseq_public_PlasmaCellsHealthy/GSE148924_raw_count_RNAseq.txt",sep="\t",header=TRUE)

#here counts are disproportionally distributed in few genes (even >2 million counts for the same gene)
#moreover, RP11-350G8.5 has only count=1 in only 2 cells. Data unusable
#scRNAseq_GSM4200472=read.table("RNAseq_public_PlasmaCellsHealthy/GSM4200472_Normal_Plasma_Cells_RawData_scRNAseq.txt",sep="\t",header=TRUE)
#scRNAseq_GSM4200472_collapsed=apply(scRNAseq_GSM4200472[,2:ncol(scRNAseq_GSM4200472)],1,sum)


scRNAseq_GSM7758185=read.table("RNAseq_public_PlasmaCellsHealthy/GSM7758185_NPCD_rep1_RSEC_ReadsPerCell_scRNAseq.csv",sep=",",header=TRUE)
scRNAseq_GSM7758185_collapsed=colSums(scRNAseq_GSM7758185[,2:ncol(scRNAseq_GSM7758185)])
scRNAseq_GSM7758186=read.table("RNAseq_public_PlasmaCellsHealthy/GSM7758186_NPCD_rep2_RSEC_ReadsPerCell_scRNAseq.csv",sep=",",header=TRUE)
scRNAseq_GSM7758186_collapsed=colSums(scRNAseq_GSM7758186[,2:ncol(scRNAseq_GSM7758186)])


scRNAseq_GSE193531=read.table("RNAseq_public_PlasmaCellsHealthy/GSE193531_umi-count-matrix_scRNAseq.csv",sep=",",header=TRUE)
nameOfGenes=scRNAseq_GSE193531$X
#subselect healthy samples (should be mix of 9 NBMs, normal bone marrow plasma cells)
scRNAseq_GSE193531_tmp=scRNAseq_GSE193531[,grepl("NBM",colnames(scRNAseq_GSE193531))]
#separate NBM samples
splittednames=sapply(strsplit(colnames(scRNAseq_GSE193531_tmp),split="\\."),"[[",4)
split_dfs = lapply(unique(splittednames), function(f) {
  df=scRNAseq_GSE193531_tmp[, splittednames == f]
  rownames(df)=nameOfGenes
  return(df)
})
names(split_dfs)=paste0("replicate_",unique(splittednames))
scRNAseq_GSE193531_collapsedList=lapply(split_dfs,rowSums)


#GSE242330_RAW_scRNAseq: cannot gunzip the replicates...
n_count_GSE148924=colSums(bulk_RNAseq_GSE148924[,(ncol(bulk_RNAseq_GSE148924)-2):ncol(bulk_RNAseq_GSE148924) ])
n_count_GSM7758185=sum(scRNAseq_GSM7758185_collapsed)
n_count_GSM7758186=sum(scRNAseq_GSM7758186_collapsed)
n_count_GSE193531=sapply(scRNAseq_GSE193531_collapsedList,sum)

bulk_RNAseq_GSE148924_filtered=bulk_RNAseq_GSE148924[bulk_RNAseq_GSE148924$ensembl==gene_id,(ncol(bulk_RNAseq_GSE148924)-2):ncol(bulk_RNAseq_GSE148924) ]
scRNAseq_GSM7758185_filtered=sum(scRNAseq_GSM7758185_collapsed["RP11.350G8.5"])
scRNAseq_GSM7758186_filtered=sum(scRNAseq_GSM7758186_collapsed["RP11.350G8.5"])
scRNAseq_GSE193531_filtered=sapply(scRNAseq_GSE193531_collapsedList,function(i){i[grepl("RP11.350G8.5",names(i))]})

#apply the formula for FPKMs (library size + gene length)
#according to https://docs.gdc.cancer.gov/Encyclopedia/pages/FPKM/:
#FPKM = [RMg * 10^9 ] / [RMt * L]
#RMg: The number of reads mapped to the gene
#RMt: The total number of read mapped to protein-coding sequences in the alignment
#L: The length of the gene in base pairs (for RP11-350G8.5 is 4236 bases ( https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000228013;r=1:154402328-154406564;t=ENST00000424435))

FPKM_GSE148924=(bulk_RNAseq_GSE148924_filtered*10^9) / (n_count_GSE148924*4236)
FPKM_GSM7758185=(scRNAseq_GSM7758185_filtered*10^9) / (n_count_GSM7758185*4236)
FPKM_GSM7758186=(scRNAseq_GSM7758186_filtered*10^9) / (n_count_GSM7758186*4236)
FPKM_GSE193531=(scRNAseq_GSE193531_filtered*10^9) / (n_count_GSE193531*4236)

df_healthy_PC=data.frame(sample_type=rep("healthy PC",14),ENSG00000228013=unlist(c(as.numeric(FPKM_GSE148924),FPKM_GSM7758185,FPKM_GSM7758186,FPKM_GSE193531)))
rownames(df_healthy_PC)=c("GSE148924_rep1","GSE148924_rep2","GSE148924_rep3","GSM7758185 scRNA","GSM7758186 scRNA",paste("GSE193531 scRNA",unique(splittednames)) )
#+ 1 and under log2
df_healthy_PC$ENSG00000228013=log2(df_healthy_PC$ENSG00000228013+1)

p2d_plotdata_details$sample_type[grepl("Baseline",p2d_plotdata_details$sample_type)]="MM patients"

colnames(p2d_plotdata_details)[2]=gene_id

p2d_plotdata_details=rbind(df_healthy_PC,p2d_plotdata_details)

print(table(p2d_plotdata_details$sample_type))
p2d <- ggboxplot(
  p2d_plotdata_details,
  x="sample_type",
  y = gene_id,
  fill = "sample_type"
) + ylab(
  "RP11-350G8.5 expression log2(FPKM+1)"
) + stat_compare_means(
  aes(x=sample_type, y = get(gene_id), group = sample_type),
  comparisons = list(
    c("MM patients","Bone Marrow"),
    c("MM patients","Normal Tissue"),
    c("MM patients","healthy PC")
  )
) + theme_bw() + theme(
  legend.position="top",
  legend.title=element_blank(),
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  legend.text=element_text(size = 16)
)






pdf("Figure_2D.pdf")
print(p2d)
dev.off()


#4*e^-10 if bone marrow+healthy PCs


#p2d_plotdata_details2=p2d_plotdata_details
#p2d_plotdata_details2$sample_type=ifelse(p2d_plotdata_details2$sample_type=="Bone Marrow","healthy PC",p2d_plotdata_details2$sample_type)





############################################################
# figure 2E
############################################################


# Add extra filter field
gene_intersections[,"is_lncRNA"] <- (
  !gene_intersections[,"is_control"] &
    !gene_intersections[,"is_ribosomal"] &
    gene_intersections[,"gene_biotype"] != "protein_coding"
)
gene_intersections[,"is_selected_by_screen_ABZB"] <- (
  gene_intersections[,"ABZB_neg.fdr"] <= 0.2
)
gene_intersections[,"is_selected_by_screen_AMO"] <- (
  gene_intersections[,"AMO_neg.fdr"] <= 0.2 
)
gene_intersections[,"is_selected_by_screen"] <- (
  gene_intersections[,"ABZB_neg.fdr"] <= 0.2 |
    gene_intersections[,"AMO_neg.fdr"] <= 0.2 
)
gene_intersections[,"is_selected_by_de"] <- (
  (gene_intersections[,"BL_NT_adj.P.Val"] <= 0.05 & gene_intersections[,"BL_NT_logFC"]>0) |
    (gene_intersections[,"BL_BM_adj.P.Val"] <= 0.05 & gene_intersections[,"BL_BM_logFC"]>0)
)


p2e_data <- gene_intersections[
# !is.na(gene_intersections[,"BL_NT_logFC"]) &
  #   !is.na(gene_intersections[,"BL_BM_logFC"]) &
  #   gene_intersections[,"BL_NT_logFC"] > 0 &
  #   gene_intersections[,"BL_BM_logFC"] > 0 &
  #   !is.na(gene_intersections[,"aggregated_score_total"]) &
  #   !gene_intersections[,"is_ribosomal"],
  # c("id", "gene_biotype", "ABZB_neg.lfc", "ABZB_neg.fdr", "AMO_neg.lfc", "AMO_neg.fdr", "aggregated_score_total")  
] %>% dplyr::arrange(dplyr::desc(aggregated_score_total))
p2e_data[,"ABZB_neg.lfc"] <- round(p2e_data[,"ABZB_neg.lfc"],2)
p2e_data[,"AMO_neg.lfc"] <- round(p2e_data[,"AMO_neg.lfc"],2)
p2e_data[,"ABZB_neg.fdr"] <- format(p2e_data[, "ABZB_neg.fdr"], digits=2,scientific = TRUE)
p2e_data[,"AMO_neg.fdr"] <- format(p2e_data[, "AMO_neg.fdr"], digits=2,scientific = TRUE)
# p2e <- grid.ftable(
#   head(p2e_data,30),
#   gp = gpar(fill = rep(c("lightgray", "gray"))),
#   each = 7
# )
p2e_data_orig=p2e_data
p2e_data[,"id"] <- ifelse(
  p2e_data[,"is_lncRNA"] &
  p2e_data[,"is_selected_by_screen"] &
  p2e_data[,"is_selected_by_de"],
  p2e_data[,"id"],
  NA
)
p2e_data <- p2e_data[!is.na(p2e_data[,"id"]),]
p2e_data <- p2e_data[order(p2e_data[,"aggregated_score_total"], decreasing = TRUE), ]
p2e_data[,"rank"] <- seq_len(nrow(p2e_data))
p2e_data <- p2e_data[p2e_data[,"rank"] <= max(p2e_data[!is.na(p2e_data[,"id"]),"rank"]),]
p2e_data[,"rank_selected"] <- ifelse(
  !is.na(p2e_data[,"id"]),
  p2e_data[,"rank"],
  NA
)

#filter the table:
#only lncRNA
#MAGECK FDR<0.2 for at least AMO1 or ABZB
#upregulated in both comparisons and significant in at least one (fold change >0 and FDR<0.05); conditions: BL_NT_logFC BL_BM_logFC
#prognostic value for OS or PFS<0.05
only_lnc_logical=p2e_data$is_lncRNA
mageck_signif_logical=as.numeric(p2e_data$ABZB_neg.fdr)<0.2 |as.numeric(p2e_data$AMO_neg.fdr) < 0.2
upreg_logical=(p2e_data$BL_NT_logFC>0 & p2e_data$BL_BM_logFC>0) & (p2e_data$BL_NT_adj.P.Val<0.05 | p2e_data$BL_BM_adj.P.Val<0.05)
prognostic_logical=p2e_data$surv_os_pval<0.05 | p2e_data$surv_pfs_pval<0.05
total_filter_logical=only_lnc_logical&mageck_signif_logical&upreg_logical&prognostic_logical

p2e_data=p2e_data[total_filter_logical,]
#p2e_data$best_adj_DGE_pval=-log10(sapply(1:length(p2e_data$BL_NT_adj.P.Val),function(i){min(p2e_data$BL_NT_adj.P.Val[i],p2e_data$BL_BM_adj.P.Val[i])}))

log10pvals_BL_NT=-log10(p2e_data$BL_NT_P.Value)
log10pvals_BL_BM=-log10(p2e_data$BL_BM_P.Value)

log10pvals_BL_NT_total=-log10(p2e_data_orig$BL_NT_P.Value)
log10pvals_BL_BM_total=-log10(p2e_data_orig$BL_BM_P.Value)
log10pvals_BL_NT_total=log10pvals_BL_NT_total[!is.na(log10pvals_BL_NT_total)&!is.infinite(log10pvals_BL_NT_total)]
log10pvals_BL_BM_total=log10pvals_BL_BM_total[!is.na(log10pvals_BL_BM_total)&!is.infinite(log10pvals_BL_BM_total)]

scaled_BL_NT=(log10pvals_BL_NT-min(log10pvals_BL_NT_total))/(max(log10pvals_BL_NT_total)-min(log10pvals_BL_NT_total))
scaled_BL_BM=(log10pvals_BL_BM-min(log10pvals_BL_BM_total))/(max(log10pvals_BL_BM_total)-min(log10pvals_BL_BM_total))

p2e_data$best_DGE_pval=sapply(1:length(scaled_BL_NT),function(i){max(scaled_BL_NT[i],scaled_BL_BM[i])})

pdf("Figure_2E.pdf")
plot(p2e_data$best_DGE_pval,p2e_data$aggregated_score_total,pch=19,xlab="Max differential expression vs. normal samples (scaled -log10 pval)",
                        ylab="Priority Score")
text(x=p2e_data$best_DGE_pval+0.01,y=p2e_data$aggregated_score_total+0.01,labels=p2e_data$gene_symbol)
#is_selected_by_screen_ABZB
#is_selected_by_screen_AMO
dev.off()

# #aggregated_score_total: is the priority score
# #min(aggregated_score_surv_os_pval/aggregated_score_surv_pfs_pval) is the max prognostic potential
# p2e <- ggplot(p2e_data) + geom_line(
#   aes(x = rank, y = aggregated_score_total),
#   linewidth = 2
# ) + geom_point(
#   aes(x = rank_selected, y = aggregated_score_total + 0.004, fill = ABZB_neg.lfc, color = "ABZB"),
#   shape = 24,
#   size = 2,
#   color = "red"
# ) + scale_fill_gradient(
#   low = "red",
#   high = "white"
# ) + new_scale_fill() + geom_point(
#   aes(x = rank_selected, y = aggregated_score_total - 0.004, fill = AMO_neg.lfc),
#   shape = 25,
#   size = 2,
#   color = "green"
# ) + scale_fill_gradient(
#   low = "darkgreen",
#   high = "white"
# ) + geom_label_repel(
#   aes(x = rank_selected, y = aggregated_score_total, label = id),
#   position = position_nudge_repel(
#     x = c(0.05, 0, -0.05, 0),
#     y = c(0.05, 0.1, -0.05, -0.1)
#   )
# ) + xlab(
#   "Rank"
# ) + ylab(
#   "Aggregated score"
# ) + theme_bw() + theme(
#   strip.text = element_text(size = 20),
#   axis.text.x = element_text(size = 20),
#   axis.text.y = element_text(size = 20),
#   axis.title.x =element_text(size = 20),
#   axis.title.y = element_text(size = 20)
# )

# pdf("Figure_2E.pdf")
# plot(p2e)
# dev.off()



############################################################
# figure 2F
############################################################


p2f <- ggplot(df_fig2F) + geom_point(
  aes(x = IL6R, y = IL6R.AS1, fill = response_status),
  size = 2,
  shape = 21
) + geom_smooth(
  aes(x = IL6R, y = IL6R.AS1, color = response_status),
) + facet_wrap(
  ~ response_status,
  nrow=3
) + theme_bw() + theme(
  legend.position = "none"
) + theme(
  strip.text = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.x =element_text(size = 20),
  axis.title.y = element_text(size = 20)
)

pdf("Figure_2F.pdf")
print(p2f)
dev.off()









