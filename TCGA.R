library("TCGAbiolinks")
library("survminer")
library("survival")
library("tidyverse")
library("DESeq2")
library("tidyverse")
library("maftools")
library("pheatmap")
library("DT")
library("SummarizedExperiment")

------------------------------------
##get a list of projects
------------------------------------
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-BRCA')

#Building a query to retrieve gene expresion data
tcga_query <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = c('TCGA-LL-A73Y-01A-11R-A33J-07','TCGA-E2-A1IU-01A-11R-A14D-07', 'TCGA-AO-A03U-01B-21R-A10J-07'))
getResults(tcga_query)


#Download data - GDC Download
GDCdownload(tcga_query)

------------------------------------------------------------------
## get gene expression data -----------
## build a query to get gene expression data for entire cohort
------------------------------------------------------------------
all_query_brca_data = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

all_query_brca_data_output <- getResults(all_query_brca_data)

#Download data - GDC Download
GDCdownload(all_query_brca_data_output)

# prepare data
tcga_query_data <- GDCprepare(all_query_brca_data, summarizedExperiment = TRUE)
  tcga_query_matrix <- assay(tcga_query_data, 'fpkm_unstrand')
tcga_query_matrix1 <- assay(tcga_query_data, 'unstranded')
 


tcga_query_matrix1[1:10,1:10]
tcga_query_matrix1[1:20,1:20]

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_query_data))
coldata <- as.data.frame(colData(tcga_query_data))

mycounts <- as.data.frame(gene_metadata)
metadata <- as.data.frame(coldata)


# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = tcga_query_matrix,
                              colData = coldata,
                              design = ~ 1)
# vst transform counts to be used in survival analysis ---------------
# Setting up countData object  
dds2 <- DESeqDataSetFromMatrix(countData = mycounts,
                              colData = coldata,
                              design = ~ 1)


dds <- DESeqDataSetFromMatrix(countData=mycounts, 
                              colData=metadata, 
                              design=~treatment, 
                              tidy=TRUE)


# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
keep1 <- rowSums(counts(dds1)) >= 10


dds <- dds1[keep1,]


# vst 
vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)
brca_matrix_vst[1:10,1:10]

idx2 <- gene_metadata %>% 
  assay %>% 
  rowVars(var = 'gene_name') %>%
  order(decreasing = TRUE) %>% 
  head(5)

idx1 <- vsd %>% 
  assay %>% 
  rowVars() %>% 
  order(decreasing = TRUE) %>% 
  head(20)

# plot
pheatmap(assay(vsd)[idx1,])
----------------------------------------------
## get 20 primary tissue sample barcodes
----------------------------------------------
  
tumor <- all_query_brca_data_output$cases[1:20]
# OR
tumor <- all_query_brca_data_output[all_query_brca_data_output$sample_type == "Primary Tumor", "cases"][1:20]
tumor
--------------------------------------------------------
# # get gene expression data from 20 primary tumors
--------------------------------------------------------
query_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumor)

# download data
GDCdownload(query_brca)

# get counts
tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, "unstranded")
brca_matrix[1:10,1:10]
brca_matrix[1:20,1:20]

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_brca_data))
coldata <- as.data.frame(colData(tcga_brca_data))

# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = brca_matrix,
                              colData = coldata,
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# vst 
vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)
brca_matrix_vst[1:10,1:10]


idx <- vsd %>% 
  assay %>% 
  rowVars(var = 'gene_id') %>%
  order(decreasing = TRUE) %>% 
  head(20)

# plot
pheatmap(assay(vsd)[idx,])


d<- duplicated(rownames(all_query_brca_data_output))
table(d)

idfound<- rownames(all_query_brca_data_output)%in% rownames(all_query_brca_data_output)
table(idfound)

head(all_query_brca_data_output)

--------------------------------------------------------------------------
## download and visualize mutation data from TCGA ----------------------
--------------------------------------------------------------------------

query_mutation <- GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Simple Nucleotide Variation',
                           access = 'open')


output_query_mutation <- getResults(query_mutation)

GDCdownload(query_mutation)

maf <- GDCprepare(query_mutation, summarizedExperiment = TRUE)
head(maf)
write_tsv(maf, "TCGA_BRCA_MAF.tsv")
-------------------------------------------------
## maftools utils to read and create dashboard
-------------------------------------------------
  
maftools.input <- read.maf(maf)
maftools.input
TCGA <- read.maf(maf = maf)


sample <- getSampleSummary(TCGA)
gene <- getGeneSummary(TCGA) 

clin <- getClinicalData(TCGA)
--------------------------------------
##Somatic interaction between genes
--------------------------------------
  
somaticinter <- somaticInteractions(maf = TCGA, top = 20, pvalue = 0.05)

driver <- oncodrive(maf = TCGA, minMut = 5, pvalMethod = "zscore")
head(driver)

---------------------------
##Protein domain changes
---------------------------
pfam <- pfamDomains(maf = TCGA, top = 10)

prSum <- pfam$proteinSummary
Dsum <- pfam$domainSummary

head(prSum)

cohort <- tcgaAvailable()
cohort <- tcgaLoad(study = "BRCA")


clin1 <- getClinicalData(cohort)

----------------------
##Survival analysiS
----------------------

clin1$days_to_last_followup[clin1$days_to_last_followup=="[Not Available]"] <- NA
clin1$vital_status[clin1$vital_status=="Alive"] <- 0
clin1$vital_status[clin1$vital_status=="Dead"] <- 1

# getting clinical data for TCGA-BRCA cohort -------------------
any(colnames(clinical) %in% c("demographic.vital_status", "diagnoses.days_to_last_follow_up", "demographic.cause_of_death"))
which(colnames(clinical) %in% c("demographic.vital_status", "diagnoses.days_to_last_follow_up", "demographic.cause_of_death"))
clinical[,c(30,65,14)]


# looking at some variables associated with survival 
table(clinical$demographic.vital_status)

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# change certain values the way they are encoded
clinical$deceased <- ifelse(clinical$demographic.vital_status == "Alive", FALSE, TRUE)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical$overall_survival <- ifelse(clinical$demographic.vital_status == "Alive",
                                    clinical$diagnoses.days_to_last_follow_up,
                                    clinical$demographic.cause_of_death)


data <- cohort@data

TCGA <- read.maf(data, clinicalData = clin1)


mafSurvival(maf = TCGA, genes = "ALS2", 
            time = "days_to_last_followup", Status = "vital_status")


mafSurvGroup(maf = TCGA, geneSet = c("","TP53"), 
             time = "days_to_last_followup", Status = "vital_status")


plotmafSummary(maf = TCGA, rmOutlier = T, addStat = "median", dashboard = T, titvRaw = F)
---------------
## oncoprint
---------------
oncoplot(maf = TCGA,
         top = 50,
         removeNonMutated = TRUE)

TCGA.titv <- titv(maf = TCGA,useSyn = TRUE)

lollipopPlot(maf = TCGA, gene = "ALS2", AACol = "HGVSp_Short",
             showMutationRate = T, showDomainLabel = T)

--------------------------------------------------------------------
## build a query to retrieve DNA methylation data --------------
--------------------------------------------------------------------
query_methly <- GDCquery(project = 'TCGA-GBM',
                         data.category = 'DNA Methylation',
                         platform = 'Illumina Human Methylation 27',
                         access = 'open',
                         data.type = 'Methylation Beta Value',
                         barcode = c('TCGA-19-0962-01B-01D-0521-05', 'TCGA-06-0137-01A-01D-0218-05'))

output_query_methyl <- getResults(query_methly)

GDCdownload(query_methly)


## plot probes showing differences in beta values between samples  

dna.meth <- GDCprepare(query_methly, summarizedExperiment = TRUE)
assay(dna.meth)  

idx <- query_brca_all_data %>% 
  assay %>% 
  rowVars() %>% 
  order(decreasing = TRUE) %>% 
  head(10)

# plot
pheatmap(assay(query_brca_all_data)[idx,])