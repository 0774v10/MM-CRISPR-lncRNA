#require objects: 
#de_data
#design.matrix




#library()

load("figure2_objects.rda")
load("figureS2_objects.rda")
load("figureS7_objects.rda")
library(survival)
library(survminer)
library(dplyr)

#read MM GDSC2
#sensitivity_MM_GDSC2=read.csv("MM_IC_Thu_Jan_25_15_41_30_2024_GDSC2.csv")
#read expressions of all cell lines (1018 unique cell lines, with COSMIC ID) https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2770/Downloads:
#download link: https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-2770/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv
expressions_total=read.table("E-MTAB-2770-query-results.tpms.tsv",sep="\t",header=TRUE)
#sensitivity pancancer for ALL possible drugs, downloaded from: https://www.cancerrxgene.org/downloads/drug_data?pathway=All
#download link: https://www.cancerrxgene.org/downloads/download/ic
sensitivity_pancancer_GDSC2_ALLdrugs=read.csv("PANCANCER_ANOVA_MonFeb5_16_04_39_2024_allGDSC2.csv")
#cohen's d and statistical tests 
df_stat_drug_sensitivity=read.table("table_exprLINC_vs_drug_sensitivity.xls",sep="\t",header=TRUE)

#subset the GDSC2, selecting only drugs involved in protein stability or degradation (including Bortezomib)
protein_stability_degradation_drugs=c("Bortezomib","CCT-018159","Elesclomol","Lenalidomide","Luminespib","MG-132","ML323","P22077","Tanespimycin")
sensitivity_pancancer_GDSC2=sensitivity_pancancer_GDSC2_ALLdrugs[sensitivity_pancancer_GDSC2_ALLdrugs$Drug.Name %in% protein_stability_degradation_drugs,]

#further subset GDSC2 for drugs in protein stability and degratation, for only MM lines
cell_lines=c("AMO-1","ARH-77","EJM","IM-9","JJN-3","KARPAS-620","KMS-11",
 "KMS-12-BM","L-363","LP-1","MM1S","MOLP-8","NCI-H929","OPM-2", 
 "RPMI-8226","SK-MM-2","U-266" )
sensitivity_MM_GDSC2=sensitivity_pancancer_GDSC2[sensitivity_pancancer_GDSC2$Cell.Line.Name %in% cell_lines,]


#define function to plot the stripchart with more features
dotplot<-function(Objectlist,col,labs,widthlines=0.2,center="median",...){
  stripchart(Objectlist,vertical=TRUE,method="jitter",xaxt="n",...)
  axis(1,at=1:length(Objectlist),labels=labs)
  for(i in 1:length(Objectlist)){
    segments(x0=i-widthlines,x1=i+widthlines,y0=median(Objectlist[[i]]),y1=median(Objectlist[[i]]),col="black")
  }
}




############################################################
# figure S7A
############################################################
gene_id <- "ENSG00000153363" #LINC00467
pS7a <- ggboxplot(do.call(rbind,lapply(
      c("bone_marrow", "normal_tissue","patient_Baseline"),
      function(sCol) cbind(
        data.frame(
          sample_type = switch(
            sCol,
            bone_marrow = "Bone Marrow",
            normal_tissue = "Normal Tissue",
            patient_Baseline = "Baseline"
          ),
          stringsAsFactors = FALSE
        ),
        t(de_data[
          gene_id,
          rownames(design.matrix)[design.matrix[,sCol]==1],
          drop = FALSE
        ])
      ))),
      x="sample_type",
      y = gene_id,
      fill = "sample_type"
    ) + ylab(
      "LINC00467 expression, log2( FPKM + 1)"
    ) + stat_compare_means(
      aes(x=sample_type, y = get(gene_id), group = sample_type),
      comparisons = list(
        c("Baseline","Bone Marrow"),
        c("Baseline","Normal Tissue")
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




pdf("Figure_S7A.pdf")
print(pS7a)
dev.off()


############################################################
# figure S7B
############################################################


response_status <- "Baseline"    
pS7b_data <- t(data_expression_patients[
  data_expression_patients[,'GENE_ID'] %in% unique(gene_intersections[,"gene_id"]) &
  data_expression_patients[,'GENE_ID'] %in% gene_list_survival,
  data_visits[data_visits[,'response_status']==response_status,'sample_id']
])
pS7b_data <- cbind(
  data_survival[
    sapply(rownames(pS7b_data),function(s) paste(strsplit(s,'_')[[1]][1:2],sep='_',collapse='_')),
    c('public_id','ttcpfs','censpfs','ttcos','censos')],
  pS7b_data
)

# Run test for all the screened genes
pS7b_coxph_all <- NULL
for (surv_type in c("os","pfs")) {    
  for (gene_id in na.omit(unique(gene_intersections[,"gene_id"]))) {
    if (gene_id %in% colnames(pS7b_data)) {
      pS7b_coxph <- suppressWarnings(coxph(formula(paste0('Surv(ttc',surv_type,', cens',surv_type,') ~ ',gene_id)), data=pS7b_data))
      pS7b_coxph_all <- rbind(pS7b_coxph_all, cbind(
        data.frame(
          gene_id = gene_id,
          surv_type = surv_type,
          coef = summary(pS7b_coxph)[["coefficients"]][,"coef"],
          exp_coef = summary(pS7b_coxph)[["coefficients"]][,"exp(coef)"],
          test_type  = c('waldtest','logtest','sctest'),
          stringsAsFactors = FALSE
        ),
        do.call(rbind,summary(pS7b_coxph)[c('waldtest','logtest','sctest')])
      ))
    }
  }
}
# write.table(
#   pS7b_coxph_all,
#   file = paste0("Figure_S1a_surfit_details.txt"),
#   row.names = FALSE,
#   col.names = TRUE,
#   sep = "\t"
# )

# Set parameters of interest
response_status <- "Baseline"    
surv_type <- "os"
gene_id <- "ENSG00000153363"

# Estimate surv for the the gene of interest
print(paste("Srrvival", surv_type))
print(do.call(rbind,summary(pS7b_coxph)[c('waldtest','logtest','sctest')]))
pS7b_data <- pS7b_data[order(pS7b_data[,gene_id], decreasing = FALSE),]

# Get best split point for Kaplan-Maier curves
for (surv_type in c("os","pfs")) {
  dPval <- data.frame()
  for (sThr in c("median","mean",1:99)) {
    if (sThr %in% c("mean","median")) {
      if (sThr=="mean") nThr <- mean(pS7b_data[,gene_id])
      if (sThr=="median") nThr <- median(pS7b_data[,gene_id])
    } else {
      nThr <- quantile(pS7b_data[,gene_id],as.integer(sThr) / 100)
    }
    nPval <- surv_pvalue(survfit(formula(paste0('Surv(ttc',surv_type,', cens',surv_type,') ~ ',gene_id,' <',nThr)) , data=pS7b_data))$pval
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
  #   file = paste0("Figure_S1a_surfit_",surv_type,".txt"),
  #   row.names = FALSE,
  #   col.names = TRUE,
  #   sep = "\t"
  # )
}

# Create surv plot
pS7b <- ggarrange(
  ggsurvplot(
    survfit(formula(paste0('Surv(ttc',"os",', cens',"os",') ~ ',gene_id,' <',25.26005)) , data=pS7b_data),
    pval = TRUE, 
    conf.int = TRUE,
    risk.table = 'abs_pct', 
    risk.table.col = "strata",
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
  # ggsurvplot(
  #   survfit(formula(paste0('Surv(ttc',"pfs",', cens',"pfs",') ~ ',gene_id,' <',6.993655)) , data=pS7b_data),
  #   pval = TRUE,
  #   conf.int = TRUE,
  #   risk.table = 'abs_pct',
  #   risk.table.col = "strata",
  #   legend.labs = c("Higer","Lower"),
  #   linetype = "strata",
  #   surv.median.line = "hv",
  #   ggtheme = theme_bw(),
  #   palette = c("#E7B800", "#2E9FDF")
  # )$plot + ylab("Survival prob. (PFS)") + theme(
  #   axis.text.x = element_text(size = 16),
  #   axis.text.y = element_text(size = 16),
  #   axis.title.x =element_text(size = 16),
  #   axis.title.y = element_text(size = 16),
  #   legend.key.size = unit(0.1, 'cm'),
  #   legend.text = element_text(size=16)
  # ),
  ncol = 1
)

pdf("Figure_S7B.pdf")
print(pS7b)
dev.off()



############################################################
# figure S7C
############################################################

#(Stripchart_IC50Bortezomib_vs_LINC00467expression_MM_splittedBy_mean_GDSC2.pdf)

############################################################
############################################################
############################################################
#1- retrieve the correct cell line name from colnames of the expression table
############################################################
############################################################
############################################################


#split colnames: name.(number). .... : first element - second element after the dot (if present)
cell_line_names_tot=sapply(colnames(expressions_total[,3:ncol(expressions_total)]),strsplit,split="\\.")
cell_line_names_firstpart=sapply(cell_line_names_tot,"[[",1)
cell_line_names_secondpart=sapply(cell_line_names_tot,"[[",2)
cell_line_names=paste(cell_line_names_firstpart,cell_line_names_secondpart,sep="-")
#remove "-" in those ending with -
dash_end=grepl("-$",cell_line_names)
tmp=cell_line_names[dash_end]
tmp2=unlist(sapply(tmp,strsplit,split="-"))
cell_line_names2=cell_line_names
cell_line_names2[dash_end]=tmp2
Xbeginning=grepl("^X",cell_line_names2)
tmp=cell_line_names2[Xbeginning]
tmp2=gsub('^.', '', tmp)
cell_line_names3=cell_line_names2
cell_line_names3[Xbeginning]=tmp2








############################################################
############################################################
############################################################
#2- Append the expression info to the sensitivity tables and filter them
############################################################
############################################################
############################################################



#extract only LINC00467 expression
LINC00467_expression=expressions_total[expressions_total$Gene.Name=="LINC00467",3:ncol(expressions_total)]

# append LINC00467 expression to table of pan cancer GDSC2
pos=match(sensitivity_pancancer_GDSC2$Cell.Line.Name,cell_line_names3)
sensitivity_pancancer_GDSC2$expression_LINC00467=rep(NA,nrow(sensitivity_pancancer_GDSC2))
sensitivity_pancancer_GDSC2$expression_LINC00467[!is.na(pos)]=as.numeric(LINC00467_expression[,pos[!is.na(pos)]])

# append LINC00467 expression to table of MM lines GDSC2
pos=match(sensitivity_MM_GDSC2$Cell.Line.Name,cell_line_names3)
sensitivity_MM_GDSC2$expression_LINC00467=rep(NA,nrow(sensitivity_MM_GDSC2))
sensitivity_MM_GDSC2$expression_LINC00467[!is.na(pos)]=as.numeric(LINC00467_expression[,pos[!is.na(pos)]])

# append LINC00467 expression to table of pan cancer GDSC2 for ALL the drugs
pos=match(sensitivity_pancancer_GDSC2_ALLdrugs$Cell.Line.Name,cell_line_names3)
sensitivity_pancancer_GDSC2_ALLdrugs$expression_LINC00467=rep(NA,nrow(sensitivity_pancancer_GDSC2_ALLdrugs))
sensitivity_pancancer_GDSC2_ALLdrugs$expression_LINC00467[!is.na(pos)]=as.numeric(LINC00467_expression[,pos[!is.na(pos)]])


#almost half of pan cancer matched with LINC expression
#almost 70% of MM lines matched with LINC expression value
#filter for only those having a non-NA LINC00467 expression
sensitivity_pancancer_GDSC2=sensitivity_pancancer_GDSC2[!is.na(sensitivity_pancancer_GDSC2$expression_LINC00467),]
sensitivity_MM_GDSC2=sensitivity_MM_GDSC2[!is.na(sensitivity_MM_GDSC2$expression_LINC00467),]
sensitivity_pancancer_GDSC2_ALLdrugs=sensitivity_pancancer_GDSC2_ALLdrugs[!is.na(sensitivity_pancancer_GDSC2_ALLdrugs$expression_LINC00467),]



############################################################
############################################################
############################################################
# 3- correlate sensitivity (IC50) to BZB with the expression value of LINC00467 and plot 
# stripchart, stratified for mean
############################################################
############################################################
############################################################
current_pancancer=sensitivity_pancancer_GDSC2
current_MM=sensitivity_MM_GDSC2
current_pancancer=current_pancancer[current_pancancer$Drug.Name %in%current_MM$Drug.Name,]

splitted_pancancer=split(current_pancancer,current_pancancer$Drug.Name)
splitted_MM=split(current_MM,current_MM$Drug.Name)

currentDrugName="Bortezomib"
sensitivity_currentdrug_pancancer=splitted_pancancer[[currentDrugName]]
sensitivity_currentdrug_MM=splitted_MM[[currentDrugName]]
cor_MM=round(cor(sensitivity_currentdrug_MM$IC50,sensitivity_currentdrug_MM$expression_LINC00467),3)
cor_pancancer=round(cor(sensitivity_currentdrug_pancancer$IC50,sensitivity_currentdrug_pancancer$expression_LINC00467),3)
tryCatch({
	cor_MM_test=round(cor.test(sensitivity_currentdrug_MM$IC50,sensitivity_currentdrug_MM$expression_LINC00467)$p.value,3)

},warning=function(x){
	cor_MM_test="not defined"
},error=function(x){cor_MM_test="not defined"})
cor_pancancer_test=round(cor.test(sensitivity_currentdrug_pancancer$IC50,sensitivity_currentdrug_pancancer$expression_LINC00467)$p.value,3)

IC50_lowexp_MM_highthresh=sensitivity_currentdrug_MM$IC50[sensitivity_currentdrug_MM$expression_LINC00467<=mean(sensitivity_currentdrug_MM$expression_LINC00467)]
IC50_highexp_MM_highthresh=sensitivity_currentdrug_MM$IC50[sensitivity_currentdrug_MM$expression_LINC00467>mean(sensitivity_currentdrug_MM$expression_LINC00467)]
#fileName=paste0("Stripchart_IC50",currentDrugName,"_vs_LINC00467expression_MM_splittedBy_mean_GDSC2.pdf")


fileName="Figure_S7C.pdf"	


pdf(fileName)
dotplot(Objectlist=list(IC50_lowexp_MM_highthresh,IC50_highexp_MM_highthresh),col=c("black","grey"),labs=c(paste("LINC00467 TPM<=",round(mean(sensitivity_currentdrug_MM$expression_LINC00467),2)),paste("LINC00467 TPM>",round(mean(sensitivity_currentdrug_MM$expression_LINC00467),2))),
	widthlines=0.2,center="median",pch=19,ylab=paste(currentDrugName,"IC50"), main=paste("sensitivity of MM cell lines to",currentDrugName))
pwilcox=round(wilcox.test(IC50_lowexp_MM_highthresh,IC50_highexp_MM_highthresh,paired=FALSE)$p.value,3)
pttest=1
tryCatch({
	pttest=round(t.test(IC50_lowexp_MM_highthresh,IC50_highexp_MM_highthresh,paired=FALSE)$p.value,3)
},warning=function(x){pttest="NA"},error=function(x){pttest="NA"})

legend("bottomright",legend=c(paste("p val Mann Withney:",pwilcox),paste("p val t test:",pttest)))
dev.off()






############################################################
# figure S7D
############################################################


pdf("Figure_S7D.pdf")
plot( df_stat_drug_sensitivity$Cohensd_standardized_effect, -log10(df_stat_drug_sensitivity$ttest_pval),pch=19,xlim=c(-1,1),ylim=c(0,0.9),col="blue",ylab="-log10 pvalue",xlab="signed effect size (Cohen's D)")
abline(v=-0.8,col="grey50")
abline(v=0.8,col="grey50")
abline(h=-log10(0.15),lty=2,col="grey20")
text(x=df_stat_drug_sensitivity$Cohensd_standardized_effect,y=-log10(df_stat_drug_sensitivity$ttest_pval)+0.03,labels=df_stat_drug_sensitivity$drug_name)

dev.off()