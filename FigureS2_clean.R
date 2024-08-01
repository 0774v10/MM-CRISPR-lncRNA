library(biomaRt)
load("figureS2_objects.rda")
load("figureS2_ABZB_counts.rda")


# lncRNA annotations downloaded from https://lncipedia.org/download
# downloaded BED files for hg19 assembly, high confidence set (protein coding genes excluded)
lncipedia=read.table("lncipedia_5_2_hc_hg19.bed")$V4

# AMO-1 RNA-Seq expression value were downloaded from https://cellmodelpassports.sanger.ac.uk/passports/SIDM00993
AMO1_expressions=read.csv("AMO1_SIDM00993_rnaseq.csv")







response_status <- "Baseline"    
p2a_data <- t(data_expression_patients[
  data_expression_patients[,'GENE_ID'] %in% unique(gene_intersections[,"gene_id"]) &
  data_expression_patients[,'GENE_ID'] %in% gene_list_survival,
  data_visits[data_visits[,'response_status']==response_status,'sample_id']
])



# mapping ENSEMBL transcript ID to ENSEMBL gene ID and gene symbol using biomaRt library
pgRNA_counts2=as.data.frame(pgRNA_counts)
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



#125 gene ids out of total 1567 IDs are not contained in patients table

#insert ENSG in the count table:
matched_ID=mapping$ensembl_gene_id[match(pgRNA_counts_filtered$ensembl_transcript_id,mapping$ensembl_transcript_id)]
matched_hgnc_symbol=mapping$hgnc_symbol[match(pgRNA_counts_filtered$ensembl_transcript_id,mapping$ensembl_transcript_id)]
pgRNA_counts_filtered$ensembl_gene_id=matched_ID
pgRNA_counts_filtered$hgnc_symbol=matched_hgnc_symbol
#filter counts without gene ID (380 out of 12102)
pgRNA_counts_filtered=pgRNA_counts_filtered[!is.na(pgRNA_counts_filtered$ensembl_gene_id),]


p2a_data_transform=log2(p2a_data+1)
#filter genes in the patients that are present in our count data
p2a_data_transform_filtered=p2a_data_transform[,colnames(p2a_data_transform)%in%pgRNA_counts_filtered$ensembl_gene_id]
#average expression across patients:
p2a_data_transform_avg=apply(p2a_data_transform_filtered,2,mean)





#create table with gene symbols and IDs for pgRNA library
splitted_norm_count=split(pgRNA_counts_filtered,pgRNA_counts_filtered$ensembl_transcript_id)
splitted_norm_count_avg=lapply(1:length(splitted_norm_count),function(i){
    current=splitted_norm_count[[i]]
    df=data.frame( ensembl_gene_id=unique(current$ensembl_gene_id),
                    symbol=unique(current$symbol),
                    hgnc_symbol=unique(current$hgnc_symbol),
                    ensembl_transcript_id=unique(current$ensembl_transcript_id))
  })

df_counts=do.call("rbind",splitted_norm_count_avg)





#transform colnames of patient data (ensID) in symbols using the count data table 
matched_patient=p2a_data_transform_avg[match(df_counts$ensembl_gene_id,names(p2a_data_transform_avg))]
df_counts$average_patient_expression=matched_patient




#set a threshold for the expression in patients
data_expression_patients_mat=as.matrix(data_expression_patients[,2:ncol(data_expression_patients)])
thresh_patients=quantile(apply(    log2(as.matrix(data_expression_patients_mat)+1)   ,1,mean,na.rm=TRUE),.25)

# #how many of the lncRNA in the wei library expressed in the MM patients?
# table(df_counts_filtered$average_patient_expression>thresh_patients)
# #97.2% (732)

# # FALSE  TRUE 
# #    21   732 
# # > thresh_patients

# #         25% 
# # 0.001319853 


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







pdf("Figure_S2A.pdf")
boxplot( log2(avg_patients_expression_allgenes+1),log2(avg_patients_expression_library+1),outline=FALSE,notch=TRUE,col=c("grey50","red"),
      ylab="Expression log2 (FPKM+1) in MM patients",names=c("all genes","wei Library"))
abline(h=thresh_patients,lty=2)
dev.off()

















#filter final table with data patients avg expression+average count
df_counts_filtered=df_counts[!is.na(df_counts$average_patient_expression),]


splittednames=sapply(lncipedia,strsplit,split="\\:")
lncRNAs_pedia=unique(unname(sapply(splittednames,"[[",1)))
#convert to ENS ID:

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






#correlate expression of AMO-1 RNA-Seq and MM patients expression
pos=match(df_counts_filtered$symbol,AMO1_expressions$symbol)
df_counts_filtered$AMO1_expression=AMO1_expressions[pos,"fpkm"]
filter=!is.na(df_counts_filtered$AMO1_expression)&!is.infinite(df_counts_filtered$AMO1_expression)&!is.na(df_counts_filtered$average_patient_expression)&!is.infinite(df_counts_filtered$average_patient_expression)

pdf("Figure_S2b.pdf")
plot(log2(df_counts_filtered$AMO1_expression+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],
      pch=19,ylab="log2 (FPKM+1) MM patients expression",xlab="log2 (FPKM+1) AMO1 expression",main=paste("wei library expression (292 shown)"))
#points(log2(df_counts_filtered$AMO1_expression+1)[df_counts_filtered$is_hit=="yes"][filter],log2(df_counts_filtered$average_patient_expression+1)[df_counts_filtered$is_hit=="yes"][filter],
      #col="red",pch=19)
corP=round(cor(log2(df_counts_filtered$AMO1_expression+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],use="complete.obs"),3)
corPtest=cor.test(log2(df_counts_filtered$AMO1_expression+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],use="complete.obs")$p.value
legend("bottomright",legend=paste0("cor Pears.: ",corP,", p=",corPtest))
dev.off()

















#use ABZB RNA-Seq count files to get RPKM and compare with MM patients as in fig S2b
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(biomaRt)
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
rownames(df_FPKM_ABZB)=sapply(strsplit(rownames(df_FPKM_ABZB),split="\\."),"[[",1)


#now join with df_counts_filtered
pos=match(df_counts_filtered$ensembl_gene_id,rownames(df_FPKM_ABZB))




df_counts_filtered$FPKM_ABZB_S1B_SCR=df_FPKM_ABZB[pos,1]
df_counts_filtered$FPKM_ABZB_S2B_SCR=df_FPKM_ABZB[pos,2]


filter=!is.na(df_counts_filtered$FPKM_ABZB_S1B_SCR)&!is.infinite(df_counts_filtered$FPKM_ABZB_S1B_SCR)&!is.na(df_counts_filtered$average_patient_expression)&!is.infinite(df_counts_filtered$average_patient_expression)
pdf("Figure_S2c.pdf")
plot(log2(df_counts_filtered$FPKM_ABZB_S1B_SCR+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],
      pch=19,ylab="log2 (FPKM+1) MM patients expression",xlab="log2 (FPKM+1) ABZB SCR 1 expression",main=paste("wei library expression (249 shown)"))
#points(log2(df_counts_filtered$FPKM_ABZB_S1B_SCR+1)[df_counts_filtered$is_hit=="yes"][filter],log2(df_counts_filtered$average_patient_expression+1)[df_counts_filtered$is_hit=="yes"][filter],
      #col="red",pch=19)
corP=round(cor(log2(df_counts_filtered$FPKM_ABZB_S1B_SCR+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],use="complete.obs"),3)
corPtest=cor.test(log2(df_counts_filtered$FPKM_ABZB_S1B_SCR+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],use="complete.obs")$p.value
legend("bottomright",legend=paste0("cor Pears.: ",corP,", p=",corPtest))
dev.off()


filter=!is.na(df_counts_filtered$FPKM_ABZB_S2B_SCR)&!is.infinite(df_counts_filtered$FPKM_ABZB_S2B_SCR)&!is.na(df_counts_filtered$average_patient_expression)&!is.infinite(df_counts_filtered$average_patient_expression)
pdf("Figure_S2d.pdf")
plot(log2(df_counts_filtered$FPKM_ABZB_S2B_SCR+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],
      pch=19,ylab="log2 (FPKM+1) MM patients expression",xlab="log2 (FPKM+1) ABZB SCR 2 expression",main=paste("wei library expression (249 shown)"))
#points(log2(df_counts_filtered$FPKM_ABZB_S2B_SCR+1)[df_counts_filtered$is_hit=="yes"][filter],log2(df_counts_filtered$average_patient_expression+1)[df_counts_filtered$is_hit=="yes"][filter],
      #col="red",pch=19)
corP=round(cor(log2(df_counts_filtered$FPKM_ABZB_S2B_SCR+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],use="complete.obs"),3)
corPtest=cor.test(log2(df_counts_filtered$FPKM_ABZB_S2B_SCR+1)[filter],log2(df_counts_filtered$average_patient_expression+1)[filter],use="complete.obs")$p.value
legend("bottomright",legend=paste0("cor Pears.: ",corP,", p=",corPtest))
dev.off()




