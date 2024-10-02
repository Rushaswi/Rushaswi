#load the TCGAbiolinks library for use
library(TCGAbiolinks)
library(SummarizedExperiment)

#get the data for all the project from GDC
GDCprojects = getGDCprojects()

#quick look at top data items
head(GDCprojects[c("project_id", "name")])

# getting the summary of data THCA(Thyroid carcinoma)
getProjectSummary("TCGA-THCA")

#specifing the data category from the summary projects 
query_TCGA = GDCquery(
  project = "TCGA-THCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts", data.type = 'Gene Expression Quantification')

#downloading the data of gdc query
GDCdownload(query = query_TCGA)

#get the cases of the downloaded samples
samplesDown <- getResults(query_TCGA,cols=c("cases"))

#  lets extract data of primary tumor
dataSmTP <- TCGAquery_SampleTypes(
  barcode = samplesDown,
  typesample = "TP"
)

# lets extract data of normal tissue
dataSmNT <- TCGAquery_SampleTypes(
  barcode = samplesDown,
  typesample = "NT"
)

# getting the query of the samples
query.samples <- GDCquery(
  project = "TCGA-THCA", 
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts", 
  barcode = c(dataSmTP, dataSmNT)
)

#downloading the query sample
GDCdownload(
  query = query.samples,
  directory = "C:/Downloads"
)

# preparing the data for the samples
dataPrep <- GDCprepare(
  query = query.samples, 
  directory = "C:/Downloads"
)

# finding the corelation above 0.6 between the samples
dataPrep_TCGA <- TCGAanalyze_Preprocessing(
  object = dataPrep, 
  cor.cut = 0.6,
  datatype = "unstranded"
) 

# normalizing the GC content from the data
library(EDASeq)
dataNorm_TCGA <- TCGAanalyze_Normalization(
  tabDF = dataPrep_TCGA,
  geneInfo = geneInfoHT,
  method = "gcContent"
)

#filetering the data unsing a method called quantile
dataFilt_TCGA <- TCGAanalyze_Filtering(
  tabDF = dataNorm_TCGA,
  method = "quantile", 
  qnt.cut =  0.25 # qnt.cut refers to quantile cutoff which 
)

# to analyse the sample data where the fdr and logfc cutoff are specified.
library(edgeR)
dataDEGs_TCGA <- TCGAanalyze_DEA(
  mat1 = dataFilt_TCGA[,dataSmNT],# read count data for one grp i.e primary tumor
  mat2 = dataFilt_TCGA[,dataSmTP],# read count data for other grp i.e normal tissue
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)  

#to save the results
write.csv(x = dataDEGs_TCGA, file = "C:/Users/dasar/OneDrive/Desktop/notes/sem2/minor project/r_script/results/sample_results2.csv")

# conversion of ensembl id to genes

library("AnnotationDbi")
library("org.Hs.eg.db")
Genelist = mapIds(org.Hs.eg.db,
                  keys=rownames(dataDEGs_TCGA), 
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")



# pathway enrichment analysis
thcaEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = Genelist
)


# to visualize the data
TCGAvisualize_EAbarplot(
  tf = rownames(thcaEA$ResBP), 
  GOBPTab = thcaEA$ResBP,
  GOCCTab = thcaEA$ResCC,
  GOMFTab = thcaEA$ResMF,
  PathTab = thcaEA$ResPat,
  nRGTab = Genelist, 
  nBar = 10,
  filename = "C:/Users/dasar/OneDrive/Desktop/notes/sem2/minor project/r_script/barplot_enrichment.pdf"
)

expressed_genes = read.csv(file = "C:/Users/dasar/OneDrive/Desktop/notes/sem2/minor project/r_script/results/sample_results.csv", stringsAsFactors = FALSE, check.names = FALSE )
head(expressed_genes)

library(ggplot2) 
str(dataDEGs_TCGA$FDR)
# Filter significantly differentially expressed genes (adjust threshold as needed)
significant_genes <- dataDEGs_TCGA[dataDEGs_TCGA$FDR < 0.05, ]


## Create volcano plot with highlighted significant genes
ggplot(significant_genes, aes(x = logFC, y = -log10(PValue), color = ifelse(FDR < 0.05 & logFC < 0, "Down-regulated", ifelse(FDR < 0.05 & logFC > 0, "Up-regulated", "Not Significant")))) +
  geom_point(size = 0.8) +
  scale_color_manual(values = c("red", "blue", "green"),
                     labels = c("Down-regulated", "Up-regulated", "Not Significant")) +  # Customize labels here
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)", title = "Volcano Plot of Differentially Expressed Genes") +
  theme_bw() +
  theme(legend.title = element_blank())  # Remove legend title



