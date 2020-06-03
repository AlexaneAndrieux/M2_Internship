library("ggplot2") 
library("FactoMineR")
library("factoextra")
library("readxl")
library("gridExtra")


#################### NORMOXIC

Data_Nx <- read_excel("~/Desktop/M2_Internship_SBRI/Reconstructions_complet.xlsx", 
                      sheet = "Nx") 


# ACP
res.pca_Nx<-PCA(Data_Nx, ncp=2) #PCA 

fviz_pca_var(res.pca_Nx, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) #Correlation circle

eig.val_Nx<-get_eigenvalue(res.pca_Nx)
fviz_eig(res.pca_Nx, addlabels = TRUE, ylim = c(0, 100)) # Eigenvalues visualization
var_Nx <- get_pca_var(res.pca_Nx) 
var_Nx$coord
var_Nx$cor
var_Nx$contrib
var_Nx$cos2

ind_Nx <- get_pca_ind(res.pca_Nx) 
ind_Nx

fviz_pca_ind(res.pca_Nx, #Graph of individuals
             geom.ind = "point", 
             palette = "Set1",
             addEllipses = TRUE, 
             legend.title = "Groups")


# Clustering 
res.hcpc_Nx <-HCPC(res.pca_Nx, method= "ward", graph=TRUE) #Clustering
plot(res.hcpc_Nx, choice="tree")#Dendrogram
Clust_Nx<-fviz_cluster(res.hcpc_Nx, #Clustering graph
                       repel = TRUE,           
                       show.clust.cent = TRUE,
                       ggtheme = theme_minimal(),
                       addEllipses=FALSE,
                       main = "Factor map") 


res.hcpc_Nx$desc.ind$para #The more typical individuals of each cluster
res.hcpc_Nx$data.clust #Individuals repartition within cluster
res.hcpc_Nx$desc.var #Cluster description according variables


