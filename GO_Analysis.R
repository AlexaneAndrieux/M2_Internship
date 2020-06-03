#check https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
#see also https://yulab-smu.github.io/clusterProfiler-book/index.html


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

BiocManager::install("org.Mm.eg.db")
# 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
BiocManager::install("pathview")
# 
BiocManager::install("RDAVIDWebService")
# 
browseVignettes("clusterProfiler")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
BiocManager::install("ReactomePA")
browseVignettes("ReactomePA")

###################################################################################################################################
# CLUSTER PROFILER ANALYSIS #######################################################################################################
###################################################################################################################################


library( clusterProfiler )
library(ggplot2) #Librairie pour le graph
library(enrichplot) #Librairie pour les graphs de ClusterProfiler
library(ReactomePA) # for Reactome Pathway Analysis
require(org.Mm.eg.db) 


#############################

setwd("~/Desktop/M2_Internship_SBRI/Analyse_enrichissement/GO_Analysis_NxHx")
data <- read.delim("Astrocytes.txt", h <- F  , sep = "\t") 
test <- as.character(EndoMarkers_Nx$...1)
ego = bitr(test, fromType <- "SYMBOL", toType <- "ENTREZID", OrgDb <- "org.Mm.eg.db") #Biological Id TRanslator
head(ego)
geneList = ego$ENTREZID
geneList

# BP analysis
ego <- enrichGO(gene          = geneList,
                universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(summary(ego))
write.csv(ego, row.names = T, file = "BPAnalysis_Astrocytes.csv")
# DotPlot BP Analysis
jpeg(filename = "graphBPAnalysis_Astrocytes.jpg", width=600, height=500)
p<-dotplot(ego, x = "qvalue", orderBy = "pvalue",
           color = "pvalue", font.size = 12, title = "", showCategory=25)
print(p)
dev.off()

#KEGG analysis
search_kegg_organism('mmu', by='kegg_code')
kk <- enrichKEGG(gene         = geneList,
                 organism     = "mmu",
                 pvalueCutoff = 0.01,
                 pAdjustMethod = "BH",
                 minGSSize = 5, 
                 maxGSSize = 500,
                 use_internal_data = FALSE)
head(kk)
write.csv(kk, row.names = T, file = "KEGGAnalysis_Astrocytes.csv")
# DotPlot KEGG Analysis
jpeg(filename = "graphKEGGAnalysis_Astrocytes.jpg", width=600, height=500)
p<-dotplot(kk, x = "qvalue", orderBy = "pvalue",
           color = "pvalue", font.size = 12, title = "", showCategory=25)
print(p)
dev.off()

#REACTOME analysis
Reac <- enrichPathway(gene = geneList,
                      organism = "mouse",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      minGSSize = 5, 
                      maxGSSize = 500,
                      readable      = TRUE)
head(Reac)
jpeg(filename = "graphREACTOMEAnalysis_Astrocytes.jpg", width=800, height=500)
p<-dotplot(Reac, x = "qvalue", orderBy = "pvalue",
           color = "pvalue", font.size = 12, title = "", showCategory=25)
print(p)
dev.off()