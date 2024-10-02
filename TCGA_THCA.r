library(TCGAbiolinks)
library(SummarizedExperiment)
GDCprojects = getGDCprojects()
head(GDCprojects[c("project_id", "name")])
getProjectSummary("TCGA-THCA")
query_TCGA = GDCquery(
  project = "TCGA-THCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts", data.type = 'Gene Expression Quantification')
GDCdownload(query = query_TCGA)
samplesDown <- getResults(query_TCGA,cols=c("cases"))
dataSmTP <- TCGAquery_SampleTypes(
  barcode = samplesDown,
  typesample = "TP"
)
dataSmNT <- TCGAquery_SampleTypes(
  barcode = samplesDown,
  typesample = "NT"
)

query.samples <- GDCquery(
  project = "TCGA-THCA", 
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts", 
  barcode = c(dataSmTP, dataSmNT)
)

GDCdownload(
  query = query.samples,
  directory = "C:/Downloads"
)

dataPrep <- GDCprepare(
  query = query.samples, 
  directory = "C:/Downloads"
)


dataPrep_TCGA <- TCGAanalyze_Preprocessing(
  object = dataPrep, 
  cor.cut = 0.6,
  datatype = "unstranded"
) 

library(EDASeq)
dataNorm_TCGA <- TCGAanalyze_Normalization(
  tabDF = dataPrep_TCGA,
  geneInfo = geneInfoHT,
  method = "gcContent"
)
dataFilt_TCGA <- TCGAanalyze_Filtering(
  tabDF = dataNorm_TCGA,
  method = "quantile", 
  qnt.cut =  0.25 # qnt.cut refers to quantile cutoff which 
)
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



