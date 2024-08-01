#require objects: 
#data_expression_patients
#data_visits
#gene_intersections
#gene_list_survival
#pgRNA_counts
#wei_library
#gstable_ABZB
#gstable_AMO

library(biomaRt)
load("figureS2_objects.rda")


# lncRNA annotations downloaded from https://lncipedia.org/download
# downloaded BED files for hg19 assembly, high confidence set (protein coding genes excluded)

lncipedia=read.table("lncipedia_5_2_hc_hg19.bed")$V4

#AMO-1 RNA-Seq expression value were downloaded from https://cellmodelpassports.sanger.ac.uk/passports/SIDM00993
AMO1_expressions=read.csv("AMO1_SIDM00993_rnaseq.csv")





nThresholdFDR=0.2
response_status <- "Baseline"    
p2a_data <- t(data_expression_patients[
  data_expression_patients[,'GENE_ID'] %in% unique(gene_intersections[,"gene_id"]) &
  data_expression_patients[,'GENE_ID'] %in% gene_list_survival,
  data_visits[data_visits[,'response_status']==response_status,'sample_id']
])


#normalize counts (lib. size and log2)
sums=apply(pgRNA_counts,2,sum,na.rm=TRUE)
normcounts=log2(t(apply(pgRNA_counts+1,1,function(i){(i/sums)*1000000})))

pgRNA_counts2=as.data.frame(normcounts)
pgRNA_counts2$symbol=wei_library$GENES
ensembl_transc_IDs=unlist(lapply(sapply(rownames(pgRNA_counts2),strsplit,split="_"),"[[",1))
#filter only for real ensembl transc_id
pgRNA_counts2$ensembl_transcript_id=ensembl_transc_IDs
pgRNA_counts_filtered=pgRNA_counts2[grepl("^ENST",pgRNA_counts2$ensembl_transcript_id),]

hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol'), 
  filters = 'ensembl_transcript_id',
  values = pgRNA_counts_filtered$ensembl_transcript_id,
  mart = hsmart
)
#125 gene ids out of 125+1442 IDs are not contained in patients table
#insert ENSG in the count table:
matched_ID=mapping$ensembl_gene_id[match(pgRNA_counts_filtered$ensembl_transcript_id,mapping$ensembl_transcript_id)]
matched_hgnc_symbol=mapping$hgnc_symbol[match(pgRNA_counts_filtered$ensembl_transcript_id,mapping$ensembl_transcript_id)]
pgRNA_counts_filtered$ensembl_gene_id=matched_ID
pgRNA_counts_filtered$hgnc_symbol=matched_hgnc_symbol
#filter counts without gene ID (380/ 380+11722)
pgRNA_counts_filtered=pgRNA_counts_filtered[!is.na(pgRNA_counts_filtered$ensembl_gene_id),]


p2a_data_transform=log2(p2a_data+1)
#filter genes in the patients that are present in our count data
p2a_data_transform_filtered=p2a_data_transform[,colnames(p2a_data_transform)%in%pgRNA_counts_filtered$ensembl_gene_id]
#average expression across patients:
p2a_data_transform_avg=apply(p2a_data_transform_filtered,2,mean)

#split and average norm counts across different sg for the same gene
splitted_norm_count=split(pgRNA_counts_filtered,pgRNA_counts_filtered$ensembl_transcript_id)
splitted_norm_count_avg=lapply(1:length(splitted_norm_count),function(i){
    current=splitted_norm_count[[i]]
    #create df with avg values and ens_ID and symbols
    mat=current[,1:8]
    avg_values=apply(mat,2,mean,na.rm=TRUE)
    df=data.frame( cbind(t(avg_values),ensembl_gene_id=unique(current$ensembl_gene_id)),
                    symbol=unique(current$symbol),
                    hgnc_symbol=unique(current$hgnc_symbol),
                    ensembl_transcript_id=unique(current$ensembl_transcript_id))
  })

df_counts=do.call("rbind",splitted_norm_count_avg)
# ENSG00000223768 ENSG00000224271 ENSG00000245750 ENSG00000253438 ENSG00000259345 
#               2               2               2               2               2 
# ENSG00000285219 ENSG00000293271 
#               2               2 

#transform colnames of patient data (ensID) in symbols using the count data table 
matched_patient=p2a_data_transform_avg[match(df_counts$ensembl_gene_id,names(p2a_data_transform_avg))]
df_counts$average_patient_expression=matched_patient
#PROBLEMS for some genes. Sometimes a single ensemblID/hgnc symbol correspond
#to multiple "original" symbols


#find the hits of the screening
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
ABZB_toplot_logical= !gstable_ABZB_rnk$is_control &!gstable_ABZB_rnk$is_ribosomal & gstable_ABZB_rnk$neg.fdr <=nThresholdFDR
AMO1_toplot_logical= !gstable_AMO_rnk$is_control &!gstable_AMO_rnk$is_ribosomal & gstable_AMO_rnk$neg.fdr <=nThresholdFDR

#significant hits to highlight in the plot
hits_ID=gstable_ABZB_rnk$ensembl_transcript_id[ABZB_toplot_logical|AMO1_toplot_logical]
pos=match(hits_ID,df_counts$ensembl_transcript_id)

#LOST 3 genes because not present in the original count table:
#hits_ID[is.na(pos)]  "ENST00000363146" (Fusion gene:GAS5-RTCA) "ENST00000437831" (LINC00624) "ENST00000607149" (RP11-337C18.8-004) are not found. 
df_counts$is_hit=rep("no",nrow(df_counts))
df_counts$is_hit[pos]="yes"

#lost further 4 hits because not found a match in the sample patient table (ENSEMBL IDs)
#df_counts[is.na(df_counts$average_patient_expression)&df_counts$is_hit=="yes",]$symbol
#"LINC00869"     "AC084809.2"    "RP11-98G7.1"   "RP11-333I13.1"

#lost 58 genes in the library that are NOT hits in our screen:
#df_counts[is.na(df_counts$average_patient_expression)&df_counts$is_hit=="no",]$symbol
#  [1] "RP1-130H16.16" "AC004893.11"   "CSAG4"         "MIR503HG"     
#  [5] "FAM211A-AS1"   "SNHG12"        "SNHG12"        "SNHG12"       
#  [9] "NEAT1"         "H19"           "FAM211A-AS1"   "CDIPT-AS1"    
# [13] "PVT1"          "RP1-122P22.2"  "TP73-AS1"      "LINC00856"    
# [17] "CTA-714B7.5"   "LRRC37A11P"    "AC007405.6"    "RP11-314C16.1"
# [21] "RP11-308D16.4" "RP11-732M18.3" "RP11-435O5.2"  "RP11-86H7.7"  
# [25] "XX-CR54.1"     "RP11-1M18.1"   "DGCR5"         "TTTY15"       
# [29] "CECR7"         "LINC01057"     "ANKRD20A19P"   "NPPA-AS1"     
# [33] "GS1-124K5.11"  "RP11-452F19.3" "NKAPP1"        "RP11-14N7.2"  
# [37] "TSTD3"         "LINC01057"     "RP11-539I5.1"  "GS1-421I3.2"  
# [41] "RP11-340E6.1"  "OR2A1-AS1"     "RP11-711M9.1"  "HULC"         
# [45] "RP11-586D19.1" "LINC00338"     "RP11-885B4.2"  "RP11-320P7.1" 
# [49] "RP11-316E14.6" "RP11-356O9.1"  "NKAPP1"        "ADAM20P1"     
# [53] "RP5-914P20.5"  "CDIPT-AS1"     "RP11-16P6.1"   "RP11-480A16.1"
# [57] "RP1-283E3.4"   "CECR7

#filter final table with data patients avg expression+average count
df_counts_filtered=df_counts[!is.na(df_counts$average_patient_expression),]
#calculate average expressions of AMO1 and ABZB across replicates
df_counts_filtered$AMO1_avg_expr=sapply(1:nrow(df_counts_filtered),function(i){mean(as.numeric(df_counts_filtered[i,3:5]),na.rm=TRUE)})
df_counts_filtered$ABZB_avg_expr=sapply(1:nrow(df_counts_filtered),function(i){mean(as.numeric(df_counts_filtered[i,6:8]),na.rm=TRUE)})



#set a threshold for the expression in both patients and cell lines
data_expression_patients_mat=as.matrix(data_expression_patients[,2:ncol(data_expression_patients)])
thresh_patients=quantile(apply(    log2(as.matrix(data_expression_patients_mat)+1)   ,1,mean,na.rm=TRUE),.25)
thresh_AMO1=quantile(df_counts_filtered$AMO1_avg_expr,.25)
thresh_ABZB=quantile(df_counts_filtered$ABZB_avg_expr,.25)
thresh_postamp=quantile(as.numeric(df_counts_filtered$PostAmpReadCount),.25)




#how many of the lncRNA in the wei library expressed in the MM patients?
table(df_counts_filtered$average_patient_expression>thresh_patients)
#97.2% (732)

# FALSE  TRUE 
#    21   732 
# > thresh_patients

#         25% 
# 0.001319853 


#Distribution plot of average expression in MM patients (all genes) vs
#the lncRNA of the library
data_expression_patients_mat=as.matrix(data_expression_patients[,2:ncol(data_expression_patients)])
avg_patients_expression_allgenes=apply(as.matrix(data_expression_patients_mat),1,mean,na.rm=TRUE)
pos=match(df_counts$ensembl_gene_id,names(avg_patients_expression_allgenes))
#(sum(is.na(pos))) -> 39 of the library not found in patients, while 776 found
avg_patients_expression_library=avg_patients_expression_allgenes[pos[!is.na(pos)]]

#remove outliers:
avg_patients_expression_allgenes=avg_patients_expression_allgenes[avg_patients_expression_allgenes<20]
avg_patients_expression_library=avg_patients_expression_library[avg_patients_expression_library<20]


# plot(density(avg_patients_expression_allgenes),lwd=2,xlab="Expression log2 (FPKM+1)",ylab="frequency")
# lines(density(avg_patients_expression_library),lwd=2,col="red")
pdf("Figure_S2A.pdf")
boxplot( log2(avg_patients_expression_allgenes+1),log2(avg_patients_expression_library+1),outline=FALSE,notch=TRUE,col=c("grey50","red"),
      ylab="Expression log2 (FPKM+1) in MM patients",names=c("all genes","wei Library"))
abline(h=thresh_patients,lty=2)
dev.off()




















#among all expressed lncRNA in patients (use 0.25 quantile expression for all genes, as before),
#how many we find in our wei library


#extract the name
splittednames=sapply(lncipedia,strsplit,split="\\:")
lncRNAs_pedia=unique(unname(sapply(splittednames,"[[",1)))
#convert to ENS ID:
library(biomaRt)
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mapping_pedia <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = lncRNAs_pedia,
  mart = hsmart
) 
pos=match(lncRNAs_pedia,mapping_pedia$hgnc_symbol)  

lncRNA_IDs_total=mapping_pedia$ensembl_gene_id[pos]
table(!is.na(lncRNA_IDs_total))
#only 3119 out of 46218+3119 lncRNA from lincpedia were correctly assigned to a ENS ID
lncRNA_IDs_total=lncRNA_IDs_total[!is.na(lncRNA_IDs_total)]
data_expression_patients_mat=as.matrix(data_expression_patients[,2:ncol(data_expression_patients)])
avg_patients_expression_allgenes=apply(as.matrix(data_expression_patients_mat),1,mean,na.rm=TRUE)

pos=match(lncRNA_IDs_total,names(avg_patients_expression_allgenes))
#lost 283 lncRNA out of 2836+283 from the total list
avg_patients_expression_lncRNAs=avg_patients_expression_allgenes[pos]
avg_patients_expression_lncRNAs=avg_patients_expression_lncRNAs[!is.na(avg_patients_expression_lncRNAs)]
thresh_patients=quantile(apply(as.matrix(data_expression_patients_mat),1,mean,na.rm=TRUE),.25)
avg_patients_expression_lncRNAs_expressed=avg_patients_expression_lncRNAs[avg_patients_expression_lncRNAs>thresh_patients]
#among all the 2836 lncRNAs in the patients (according to lincipedia and biomart),
#2238 are "expressed" (>25% quantile)
pos=match(names(avg_patients_expression_lncRNAs_expressed),gene_intersections$gene_id[gene_intersections$gene_biotype!="protein_coding"])

#221 of the patients lncRNA expressed are found in wei library(excluding protein_coding)
#2017 of the lncRNA expressed in patients are not in our wei library







#correlate expression of AMO-1 RNA-Seq and MM patients expression

pos=match(df_counts_filtered$symbol,AMO1_expressions$symbol)

df_counts_filtered$AMO1_expression=AMO1_expressions[pos,"fpkm"]

#only 1/third of the library used was found in the expression table
# FALSE  TRUE 
#   292   461 
# filter=!is.na(log2(df_counts_filtered$AMO1_expression))&!is.infinite(log2(df_counts_filtered$AMO1_expression))&!is.na(log2(df_counts_filtered$average_patient_expression))&!is.infinite(log2(df_counts_filtered$average_patient_expression))
# pdf("extra_figures/Scatterplot_expr_AMO1rnaseq_vs_MMpatients.pdf")
# plot(log2(df_counts_filtered$AMO1_expression)[filter],log2(df_counts_filtered$average_patient_expression)[filter],
#       pch=19,ylab="log2 FPKM MM patients expression",xlab="log2 FPKM AMO1 expression",main=paste("wei library expression (246 shown)"))
# points(log2(df_counts_filtered$AMO1_expression)[df_counts_filtered$is_hit=="yes"][filter],log2(df_counts_filtered$average_patient_expression)[df_counts_filtered$is_hit=="yes"][filter],
#       col="red",pch=19)
# corP=round(cor(log2(df_counts_filtered$AMO1_expression)[filter],log2(df_counts_filtered$average_patient_expression)[filter],use="complete.obs"),3)
# corPtest=cor.test(log2(df_counts_filtered$AMO1_expression)[filter],log2(df_counts_filtered$average_patient_expression)[filter],use="complete.obs")$p.value
# legend("topleft",legend=c(paste0("NOT hit, cor:",corP,"; p=",round(corPtest,3)),"hit"),col=c("black","red"),pch=19)
# dev.off()


filter=!is.na(df_counts_filtered$AMO1_expression)&!is.infinite(df_counts_filtered$AMO1_expression)&!is.na(df_counts_filtered$average_patient_expression)&!is.infinite(df_counts_filtered$average_patient_expression)
pdf("Figure_S2b.pdf")
plot(log2(df_counts_filtered$AMO1_expression+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],
      pch=19,ylab="log2 (FPKM+1) MM patients expression",xlab="log2 (FPKM+1) AMO1 expression",main=paste("wei library expression (292 shown)"))
#points(log2(df_counts_filtered$AMO1_expression+1)[df_counts_filtered$is_hit=="yes"][filter],log2(df_counts_filtered$average_patient_expression+1)[df_counts_filtered$is_hit=="yes"][filter],
      #col="red",pch=19)
corP=round(cor(log2(df_counts_filtered$AMO1_expression+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],use="complete.obs"),3)
corPtest=cor.test(log2(df_counts_filtered$AMO1_expression+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],use="complete.obs")$p.value
dev.off()






#use ABZB RNA-Seq count files to get RPKM and compare with MM patients as in fig S2b
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(biomaRt)
load("figureS2_ABZB_counts.rda")
Genes=genes(TxDb.Hsapiens.UCSC.hg19.knownGene)


mart <- useMart("ENSEMBL_MART_ENSEMBL")#,host = "asia.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "entrezgene_id",
    "ensembl_gene_id",
    "external_gene_name"),
  filter = "entrezgene_id",
  values = Genes$gene_id, uniqueRows=TRUE)

pos=match(Genes$gene_id,annotLookup$entrezgene_id)


Genes$ensembl_id=annotLookup$ensembl_gene_id[pos]



#match with counts 

ensembl_id_counts=sapply(strsplit(rownames(rnaseq_data_dgelist$counts),split="\\."),"[[",1)

pos=match(ensembl_id_counts,Genes$ensembl_id)

#genes_final=width(Genes)[pos]
counts_ABZB_filtered=rnaseq_data_dgelist$counts[!is.na(pos),]
Genes_filtered=Genes[pos[!is.na(pos)]]


#select only SCR ABZB samples (positions 3 and 6 for SCR replicates)
counts_ABZB_filtered=counts_ABZB_filtered[,c(3,6)]

norm_factors_filtered=calcNormFactors(rnaseq_data_dgelist)$samples$lib.size[c(3,6)]

#apply the formula for FPKMs (library size + gene length)
#according to https://docs.gdc.cancer.gov/Encyclopedia/pages/FPKM/:
#FPKM = [RMg * 109 ] / [RMt * L]
#RMg: The number of reads mapped to the gene
#RMt: The total number of read mapped to protein-coding sequences in the alignment
#L: The length of the gene in base pairs

FPKM_1=(counts_ABZB_filtered[,1]*10^9) / (norm_factors_filtered[1]*width(Genes_filtered))
FPKM_2=(counts_ABZB_filtered[,2]*10^9) / (norm_factors_filtered[2]*width(Genes_filtered))


df_FPKM_ABZB=cbind(FPKM_1,FPKM_2)
colnames(df_FPKM_ABZB)=colnames(counts_ABZB_filtered)







#export raw count files for ABZB SCR replicates and KO replicates
pos_abzb=grepl("^S.B",colnames(rnaseq_data_dgelist$counts))
ABZB_counts=rnaseq_data_dgelist$counts[,pos_abzb]

write.table(ABZB_counts,file="ABZB_raw_counts.tsv",sep="\t",quote=FALSE,row.names=TRUE,col.names=NA)

