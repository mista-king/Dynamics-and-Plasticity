library(TCGAbiolinks)
TCGAbiolinks:::getGDCprojects()$project_id
cancer_type="TCGA-LUAD"
clinical=GDCquery_clinic(project=cancer_type,type="clinical")

dim(clinical)
clinical[1:4,1:4]
head(colnames(clinical))
query <- GDCquery(project = cancer_type, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query, files.per.chunk = 100)
expdat <- GDCprepare(query = query)

library(SummarizedExperiment)
tpm <- expdat@assays@data$tpm_unstrand
dimnames(tpm) <- list(expdat@rowRanges$gene_name,expdat@colData@rownames)


meta = clinical
#同样以第一列作为列名
meta <- column_to_rownames(meta,var = "submitter_id")
colnames(meta)
#筛选我们感兴趣的临床信息
meta=meta[,colnames(meta) %in% c("vital_status",
                                 "days_to_last_follow_up",
                                 "days_to_death",
                                 "race",
                                 "gender",
                                 "age_at_index",
                                 "tumor_stage")]
dim(meta) #585个
head(colnames(tpm))
head(rownames(meta))
#调整、筛选临床样本信息数量与顺序与表达矩阵相同
meta=meta[match(substr(colnames(tpm),1,12),rownames(meta)),]
dim(tpm)
head(colnames(tpm))
dim(meta)
head(rownames(meta))


#1、计算生存时间
meta$days_to_death[is.na(meta$days_to_death)] <- 0   #缺失值标记为0
meta$days_to_last_follow_up[is.na(meta$days_to_last_follow_up)] <- 0
meta$days=as.numeric(meta[,1])+as.numeric(meta[,6])
#时间以月份记，保留两位小数
meta$time=round(meta$days/30,2)
#2、根据生死定义活着是0，死的是1
meta$event=ifelse(meta$vital_status=='Alive',0,1)
table(meta$event)
#3 年龄分组(部分样本缺失，考虑可能的影响应该不大)
meta$age_at_index[is.na(meta$age_at_index)] <- 0
meta$age_at_index=as.numeric(meta$age_at_index)
meta$age_group=ifelse(meta$age_at_index>median(meta$age_at_index),'older','younger')
table(meta$age_group)
#4 癌症阶段
table(meta$tumor_stage)
#5 race 人种
table(meta$race)
#6 性别 gender
table(meta$gender)




library(survival)
library(survminer)
exprSet=tpm  #套流程
g = rownames(exprSet)[c(7197,6556,776)]# 随便选一个
meta$gene = ifelse(exprSet[g,]>median(exprSet[g,]),'high','low')
sfit1=survfit(Surv(time, event)~gene, data=meta)
ggsurvplot(sfit1,pval =TRUE, data = meta, risk.table = TRUE)




















library(TCGAbiolinks)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(survival)
library(survminer)

TCGAbiolinks:::getGDCprojects()$project_id
cancer_type="TCGA-LUAD"
clinical=GDCquery_clinic(project=cancer_type,type="clinical")

dim(clinical)
clinical[1:4,1:4]
head(colnames(clinical))
query <- GDCquery(project = cancer_type, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
GDCdownload(query, files.per.chunk = 100)
expdat <- GDCprepare(query = query)
count_matrix <- assay(expdat)
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(count_matrix),typesample = c("NT"))
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(count_matrix),typesample = c("TP"))
rownames(count_matrix)=gsub("\\.\\d*", "", rownames(count_matrix))

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(count_matrix)
hsa_symbols<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',"chromosome_name", "start_position","end_position", "band"),filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)

gexp <- count_matrix[c('ENSG00000141510'),samplesTP]
names(gexp) <-  sapply(strsplit(names(gexp),'-'),function(x) paste0(x[1:3],collapse="-"))
clinical$'ENSG00000141510' <- gexp[clinical$submitter_id]
df<-subset(clinical,select =c(submitter_id,vital_status,days_to_death,ENSG00000141510))
df <- df[!is.na(df$ENSG00000141510),]
df$exp <- ''
df[df$ENSG00000141510 >= mean(df$ENSG00000141510),]$exp <- "H"
df[df$ENSG00000141510 <  mean(df$ENSG00000141510),]$exp <- "L"
df[df$vital_status=='Dead',]$vital_status <- 2
df[df$vital_status=='Alive',]$vital_status <- 1
df$vital_status <- as.numeric(df$vital_status)
fit <- survfit(Surv(days_to_death, vital_status)~exp, data=df)
surv_pvalue(fit)$pval.txt 
ggsurvplot(fit,pval=TRUE)











library(cgdsr)
library(DT)
mycgds <- CGDS("http://www.cbioportal.org/")
test(mycgds) 
all_TCGA_studies <- getCancerStudies(mycgds)
DT::datatable(all_TCGA_studies)
blca2018='blca_tcga_pan_can_atlas_2018'
luad2018='luad_tcga_pan_can_atlas_2018'
lihc2018='lihc_tcga_pan_can_atlas_2018'
skcm2018='skcm_tcga_pan_can_atlas_2018'
prad2018='prad_tcga_pan_can_atlas_2018'
chol2018='chol_tcga_pan_can_atlas_2018'
thca2018='thca_tcga_pan_can_atlas_2018'
esca2018='esca_tcga_pan_can_atlas_2018'
brca2018='brca_tcga_pan_can_atlas_2018'
paad2018='paad_tcga_pan_can_atlas_2018'
stad2018='stad_tcga_pan_can_atlas_2018'
kirc2018='kirc_tcga_pan_can_atlas_2018'
coad2018='coadread_tcga_pan_can_atlas_2018'
ov2018='ov_tcga_pan_can_atlas_2018'
getCaseLists(mycgds,brca2018)[,c(1,2)]
getGeneticProfiles(mycgds,brca2018)[,c(1,2)]



#B_Macrophage
choose_genes1 <- c("TGFB1","AGT",'CCL3','IFNG','TNF')   #index1
choose_genes2 <- c("CSF1","PTPRC",'SELP','CD40LG','IL15','ADAM17','CXCL12','APOE','APP','HAS2','TNFSF13B','LAMA2','IGF1','IL6')   #index2
choose_genes3 <- c("ICAM4","LAMB1",'HSPG2','AGRN','THBS1','A2M','FBN1','GMFB','PCDH1')   #index3
choose_genes4 <- c("CCL4","TGFB3",'THBS2','CD34','CCL28','HMGB1')   #index4
choose_genes5 <- c("MIF","CCL15",'CELSR1','MDK','PTN','VEGFA','CADM1','JAG1')   #index5
choose_genes6 <- c("COL4A1","ITGAM",'TNFSF13','ICAM1','EDN1','ANXA1','SPP1','CXCL2','FGF2','SEMA4D','XCL2','CEACAM1','SAA1','EDN2','NECTIN1')   #index6
choose_genes7 <- c("BMP7","HLA-DRA",'COL1A1','CX3CL1','PIK3CB','CCL20','CXCL14','HBEGF','HMGB2','ALOX5AP','CXCL3','NMU','CXCL13','NPW','F11R','CXCL16','HEBP1')   #index7
choose_genes8 <- c("CCL2","CD28",'HLA-DMA','CLCF1','IL34','PDGFB','SPN')   #index8

#NK_CD4T
choose_genes1 <- c("TGFB1","JAG1",'ADAM17','TNFSF13B')   #index1
choose_genes2 <- c("PTPRC","IFNG",'IL15','HAS2','CXCL12','APP','TNF','HMGB2','MIF','CALR','IGF1')   #index2
choose_genes3 <- c("YARS","CELSR1",'NECTIN1','LAMA2','FGF2')   #index3
choose_genes4 <- c("OSM","CCL3",'LIF','CCL4','CCL2')   #index4
choose_genes5 <- c("EDN1","CD40",'TNFSF13','ITGB1','CXCL2','CXCL16','CCL5','LAMB2','HLA-E','PROS1','JAM2')   #index5
choose_genes6 <- c("HMGB1","CD28",'SPP1')   #index6
choose_genes7 <- c("TNFSF12","PTDSS1",'HLA-F','PLAU','HLA-DRA','SELE','ANXA1','NAMPT','ARF1','PTPRF','CLSTN1','SEMA4D','DUSP18')   #index7
choose_genes8 <- c("CALM1","TNC",'NPNT','A2M','CADM3','SECTM1','PTN','VCAM1','THBS1')   #index8
choose_genes9 <- c("GSTP1","CD2",'XCL2','ITGB7','ICAM3','COL1A1','HSPG2','THBS2','NMU','LAMB1','CLDN11')   #index9
choose_genes10 <- c("DSG2","CDH1",'COPA','FBN1','CX3CL1','FBN1','COL4A1','FN1','COL4A1','CD58','PCDH1','GMFB')   #index10


#CD4
choose_genes1 <- c("COL4A1","LAMC2",'CELSR1','LIF','GMFB','JAM3','THBS1','CADM3','NECTIN3','PCDH1','GRN','TGFB3')   #index1
choose_genes2 <- c("C3","GSTP1",'HLA-A','ITGA4','VCAM1','LTB','NECTIN1')   #index2
choose_genes3 <- c("IGF1","CCL20",'PECAM1','TNFSF13')   #index3
choose_genes4 <- c("CCL11","CLSTN1",'IL1RN','DUSP18','IL15','IL16','PTPRC','AGT','CD28','IL18','ICAM2','PIK3CB','SELL','CXCL10','IL33','RSPO3')   #index4
choose_genes5 <- c("NMU","TNFSF12",'CXCL13','ICAM1','CXCL14','HAS2')   #index5
choose_genes6 <- c("LGALS3","SAA1",'HMGB1','TNF','CD2','CEACAM1','TGFB1','ALCAM','XCL2','CXCL3','CALM1','CALM1','TNC','CXCL16','NMB')   #index6


#Macrophage_Neutrophil
choose_genes1 <- c("FGF2","VEGFA",'TGFB1','IFNG','TNF','AGT','APOE','CSF1','CCL3','HAS2','IL1B','CCL2','ICAM1')   #index1
choose_genes2 <- c("VWF","CCL28",'PGF','NMB','XCL2','CXCL16','NPW','ITGAM','SAA1')   #index2
choose_genes3 <- c("EDN1","SPP1",'ADAM17','CXCL12','IL33','HMGB1','NAMPT')   #index3
choose_genes4 <- c("ANGPT2","CRLF1",'HLA-DMA','CLCF1','COL4A1')   #index4
choose_genes5 <- c("THBS2","CD40LG",'IL15','IGF1','APP','LAMA2','IL6','OSM')   #index5
choose_genes6 <- c("CALR","COL5A3",'INHBA','OCLN','LAMB2','JAG1','MIF','NPNT')   #index6
choose_genes7 <- c("SERPING1","IL34",'MDK','LIF','AREG','VEGFC','CEACAM1','CCL15','SEMA3C','PIK3CB','CCL5','PDGFB')   #index7
choose_genes8 <- c("WNT5A","GMFB",'CXCL14','CD58','TFPI','FN1','HBEGF')   #index8





#mRNA
mycaselist <- getCaseLists(mycgds,blca2018)[6,1] 
mygeneticprofile <- getGeneticProfiles(mycgds,blca2018)[9,1]
expr <- getProfileData(mycgds,choose_genes1,mygeneticprofile,mycaselist)

#mutation
mycaselist <- getCaseLists(mycgds,ov2018)[8,1] 
mygeneticprofile <- getGeneticProfiles(mycgds,ov2018)[6,1]
mut_df <- getProfileData(mycgds,choose_genes,mygeneticprofile,mycaselist)

mut_df <- apply(mut_df,2,as.factor)
#去除“NAN”
mut_df[mut_df == "NaN"] = ""
#去除“NA”
mut_df[is.na(mut_df)] = ""
mut_df[mut_df != ''] = "MUT"


#cna
mycaselist <- getCaseLists(mycgds,ov2018)[3,1] 
mygeneticprofile <- getGeneticProfiles(mycgds,ov2018)[3,1]
cna <- getProfileData(mycgds,choose_genes,mygeneticprofile,mycaselist)

#clinical
mycaselist <- getCaseLists(mycgds,blca2018)[1,1] 
myclinicaldata <- getClinicalData(mycgds,mycaselist)
choose_columns=c('AGE','DFS_MONTHS','DFS_STATUS',
                 'NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT',
                 'OS_STATUS','OS_MONTHS','RADIATION_THERAPY','SEX'
)
choose_clinicaldata <- myclinicaldata[,choose_columns]

library(survival)
dat <- choose_clinicaldata[choose_clinicaldata$OS_MONTHS >0,]
table(dat$OS_STATUS)
my.surv <- Surv(dat$OS_MONTHS,dat$OS_STATUS=='1:DECEASED') 
kmfit1 <- survfit(my.surv~1) 
library("survminer")
ggsurvplot(kmfit1,data = dat)



length(intersect(rownames(choose_clinicaldata),rownames(mut_df)))
dat2 <- cbind(choose_clinicaldata[rownames(mut_df),c('OS_STATUS','OS_MONTHS')],mut_df)
dat2 <- dat2[dat2$OS_MONTHS > 0,]
dat2 <- dat2[!is.na(dat$OS_STATUS),]
dat2$OS_STATUS <- as.character(dat2$OS_STATUS) 
attach(dat2)
my.surv <- Surv(OS_MONTHS,OS_STATUS=='1:DECEASED')
kmfit4 <- survfit(my.surv~TP53,data = dat2)  
kmfit5 <- survfit(my.surv~RPL13A,data = dat2)  
ggsurvplot(kmfit4,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE,
             title="Survival plot by TP53 Mutaiton")
ggsurvplot(kmfit5,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE,
             title="Survival plot by RPL13A Mutaiton") 
detach(dat2)



dat3=cbind(choose_clinicaldata[,c('OS_STATUS','OS_MONTHS')],expr[rownames(choose_clinicaldata),])
dat3=dat3[dat3$OS_MONTHS > 0,]
dat3=dat3[!is.na(dat3$OS_STATUS),]
dat3$OS_STATUS=as.character(dat3$OS_STATUS)
dat3$total=rowSums(dat3[,c(3:7)])

p <- ggboxplot(dat3, x="OS_STATUS", y="total", color = "OS_STATUS",
                 palette = "jco", add = "jitter")
p+stat_compare_means(method = "t.test")
dat3=dat3[!is.na(dat3$total),]
dat3$total = ifelse(dat3$total > median(dat3$total),'high','low')
attach(dat3)
table(dat3$total)
my.surv <- Surv(dat3$OS_MONTHS,dat3$OS_STATUS=='1:DECEASED')
kmfit1 <- survfit(my.surv~total,data = dat3) 
ggsurvplot(kmfit1,conf.int =F, pval = T)


getCaseLists(mycgds,m)[,c(1,2)]
getGeneticProfiles(mycgds,m)[,c(1,2)]
pdf('THCA_MacroNeu.pdf')
for (i in c(1:8) ) {
  j=paste('choose_genes',i,sep='')
  m=ov2018
  mycaselist <- getCaseLists(mycgds,m)[6,1] 
  mygeneticprofile <- getGeneticProfiles(mycgds,m)[12,1]
  expr <- getProfileData(mycgds,get(j),mygeneticprofile,mycaselist)
  mycaselist <- getCaseLists(mycgds,m)[1,1] 
  myclinicaldata <- getClinicalData(mycgds,mycaselist)
  choose_columns=c('OS_STATUS','OS_MONTHS','SEX'
  )
  choose_clinicaldata <- myclinicaldata[,choose_columns]
  dat3=cbind(choose_clinicaldata[,c('OS_STATUS','OS_MONTHS')],expr[rownames(choose_clinicaldata),])
  dat3=dat3[dat3$OS_MONTHS > 0,]
  dat3=dat3[!is.na(dat3$OS_STATUS),]
  dat3$OS_STATUS=as.character(dat3$OS_STATUS)
  dat3$total=rowSums(dat3[,c(3:ncol(dat3))])
  dat3=dat3[!is.na(dat3$total),]
  dat3$total = ifelse(dat3$total > median(dat3$total),'high','low')
  attach(dat3)
  table(dat3$total)
  my.surv <- Surv(dat3$OS_MONTHS,dat3$OS_STATUS=='1:DECEASED')
  kmfit1 <- survfit(my.surv~total,data = dat3) 
  p=ggsurvplot(kmfit1,conf.int =F, pval = T)
  print(p)
}
dev.off()
  