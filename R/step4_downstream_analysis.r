library(Seurat)
library(dplyr)
library(ggsci)
library(tidyverse)
library(patchwork)
library(scPlant)
library(monocle3)

rm(list=ls())

########################################***Monocle Analysis***########################################
seu_obj <- readRDS("./R_result_renew.rds")
seu_ref <- readRDS("../data/SeuratObj.rds")

########################################***Cell type anno***########################################
seu_ref <- readRDS(ref_file)
seu_ref$label <- seu_ref$celltype
sce_ref <- Seurat::as.SingleCellExperiment(seu_ref)  # Must transforme to sce obj
seu_obj <- AutoAnnotate_SingleR(seu_obj, sce_ref, ref_type="single-cell")  #cell type auto anno by SingleR

plot1<-DimPlot(seu_obj, reduction = "umap", group.by='predicted_label') + scale_color_npg()
plot2<-DimPlot(seu_obj, reduction = "umap", label = TRUE, repel=TRUE) + scale_color_npg()
plot3<-DimPlot(seu_obj, reduction = "umap", split.by='orig.ident', group.by='predicted_label', repel=TRUE) + scale_color_npg()
ggsave("umap_anno.pdf", plot=plot1, width = 8, height = 6)
ggsave("umap_pre.pdf", plot=plot2, width = 8, height = 6)
ggsave("umap_split.pdf", plot=plot3, width = 14, height = 4)



########################################***Trajectory analysis***########################################

#####################################***data trasn to cds***#####################################
data<-GetAssayData(seu_obj,assay ='RNA',slot ='counts')
cell_metadata <-seu_obj@meta.data
gene_annotation <-data.frame(gene_short_name =rownames(data))
rownames(gene_annotation)<-rownames(data)

cds <-new_cell_data_set(data,
          cell_metadata =cell_metadata,
          gene_metadata =gene_annotation)


#####################################***data preprocess***#####################################

cds <- preprocess_cds(cds, num_dim=50)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, preprocess_method="PCA")
plot_cells(cds)

#从Seurat对象导入UMAP坐标，并通过monocle3绘图函数进行绘制
cds.embed<- cds@int_colData$reducedDims$UMAP   #获取CDS对象中umap的坐标
int.embed <- Embeddings(seu_obj, reduction = "umap")   #获取Seurat对象中umap的坐标
###########这里之所以不可以拿Seurat对象中的UMAP是因为细胞ID的数量不同，毕竟是两种降维聚类方法，所以我们要提取坐标中细胞ID一致的，然后绘图，才可以更好比较结果
int.embed <- int.embed[rownames(cds.embed),]#取交集
cds@int_colData$reducedDims$UMAP <- int.embed#赋值

p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters")


#####################################***data analysis***#####################################

cds <- cluster_cells(cds, reduction_method="UMAP", cluster_method="louvain")
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "predicted_label",
               label_groups_by_cluster = FALSE, 
               label_leaves = FALSE, 
               label_branch_points = FALSE)
               
p5 <- plot_cells(cds, color_cells_by = "predicted_label",
                 label_groups_by_cluster=FALSE,
                 label_leaves=TRUE,
                 label_branch_points=TRUE,
                 group_label_size=1.5)
                 
cds <- order_cells(cds)  
## options(browser="firefox")  如果弹不出交互窗口，更改浏览器即可。

p6 <- plot_cells(cds,
          color_cells_by = "pseudotime",
          label_cell_groups=FALSE,
          label_leaves=FALSE,
          label_branch_points=FALSE,
          graph_label_size=1.5)

subset_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)

top10 <- subset_pr_test_res %>% top_n(n=10, morans_I) %>%
 pull(gene_short_name) %>% as.character()
p7 <- plot_genes_in_pseudotime(cds[top10,], color_cells_by="predicted_label", 
                             min_expr=0.5, ncol = 2)



########################################***WGCNA Analysis***########################################

DO_WGCNA <- function(seu_obj, type="unsigned", cortype="pearson"){
    library(WGCNA)
    library(Seurat)
    library(tidyverse)
    library(reshape2)
    library(stringr)

    datadf <- as.matrix(seu_obj@assays$RNA@data )
    idd1 <- seu_obj@meta.data
    Inter.id1<-cbind(rownames(idd1),idd1$seurat_clusters)
    rownames(Inter.id1)<-rownames(idd1)
    colnames(Inter.id1)<-c("CellID","Celltype") ## Need Change 
    Inter.id1<-as.data.frame(Inter.id1)
    head(Inter.id1)
    Inter1<-datadf[,Inter.id1$CellID]
    Inter2<-as.matrix(Inter1)
    Inter2[1:4,1:4]

    # Create the pseudocells 
    pseudocell.size = 10 ## 10 test
    new_ids_list1 = list()
    length(levels(factor(Inter.id1$Celltype)))

    for (i in 1:length(levels(factor(Inter.id1$Celltype)))) {
        cluster_id = levels(Inter.id1$Celltype)[i]
        cluster_cells <- rownames(Inter.id1[Inter.id1$Celltype == cluster_id,])
        cluster_size <- length(cluster_cells)     
        pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
        pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
        names(pseudo_ids) <- sample(cluster_cells)    
        new_ids_list1[[i]] <- pseudo_ids      
    }

    new_ids <- unlist(new_ids_list1)
    new_ids <- as.data.frame(new_ids)
    head(new_ids)
    new_ids_length <- table(new_ids)
    new_ids_length

    new_colnames <- rownames(new_ids)  ###add
    #rm(all.data1)
    gc()
    colnames(datadf)  
    all.data<-datadf[,as.character(new_colnames)] ###add
    all.data <- t(all.data)###add
    new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                        list(name=new_ids[,1]),FUN=mean)
    rownames(new.data)<-new.data$name
    new.data <- new.data[,-1]
    new_ids_length<-as.matrix(new_ids_length)##

    short<-which(new_ids_length<10)##
    new_good_ids<-as.matrix(new_ids_length[-short,])##
    result<-t(new.data)[,rownames(new_good_ids)]
    dim(result)

    seu_obj <- FindVariableFeatures(seu_obj,nfeatures = 5000)
    colnames(result)[grepl("[12]_Cel",colnames(result))]
    Cluster1 <- result[intersect(Seurat::VariableFeatures(seu_obj),rownames(result)),]


    type <- "unsigned"  # 官方推荐 "signed" 或 "signed hybrid"
    corType <- "pearson" # 相关性计算  官方推荐 biweight mid-correlation & bicor  corType: pearson or bicor 
    corFnc <- ifelse(corType=="pearson", cor, bicor)
    # corFnc
    maxPOutliers <- ifelse(corType=="pearson",1,0.05) # 对二元变量，如样本性状信息计算相关性时，或基因表达严重依赖于疾病状态时，需设置下面参数关联样品性状的二元变量时，设置
    robustY <- ifelse(corType=="pearson",T,F)
    dataExpr <- as.matrix(Cluster1)

    m.mad <- apply(dataExpr,1,mad)
    dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

    ## 转换为样品在行，基因在列的矩阵
    dataExpr <- as.data.frame(t(dataExprVar))
    dim(dataExpr)
    head(dataExpr)[,1:8]

    ## 检测缺失值
    gsg = goodSamplesGenes(dataExpr, verbose = 3)
    # gsg$allOK
    # gsg$goodSamples

    if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
        printFlush(paste("Removing genes:", 
                        paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", 
                        paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
    }

    nGenes = ncol(dataExpr)
    nSamples = nrow(dataExpr)

    ## 查看是否有离群样品
    pdf("sampleTree.pdf")
    sampleTree = hclust(dist(dataExpr), method = "average")
    plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
    dev.off()

    powers <- c(c(1:10), seq(from = 12, to=30, by=2))
    sft <- pickSoftThreshold(dataExpr, powerVector=powers, 
                            networkType="signed", verbose=5)
    
    powers <- sft$powerEstimate

    cor <- WGCNA::cor

    net <- blockwiseModules(dataExpr, power = powers, maxBlockSize = nGenes,
                           TOMType = "unsigned", minModuleSize = 10,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs=TRUE, corType = corType, 
                           maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                           saveTOMFileBase = paste0("dataExpr", ".tom"),
                           verbose = 3)
    
    moduleLabels <- net$colors
    moduleColors <- labels2colors(moduleLabels)
    # moduleColors
    # Plot the dendrogram and the module colors underneath
    # 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    
    # module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
    MEs = net$MEs

    ### 不需要重新计算，改下列名字就好
    ### 官方教程是重新计算的，起始可以不用这么麻烦
    MEs_col = MEs
    colnames(MEs_col) = paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
    MEs_col = orderMEs(MEs_col)

    # 根据基因间表达量进行聚类所得到的各模块间的相关性图
    # marDendro/marHeatmap 设置下、左、上、右的边距
    head(MEs_col)
    plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                          marDendro = c(3,3,2,4),
                          marHeatmap = c(3,4,2,2),
                          plotDendrograms = T,
                          xLabelsAngle = 90)
    

    which.module="turquoise"; 
    ME=mergedMEs[, paste("ME",which.module, sep="")]
    par(mfrow=c(2,1), mar=c(0,4.1,4,2.05))
    plotMat(t(scale(dataExpr[,moduleColors==which.module ]) ),
            nrgcols=30,rlabels=F,rcols=which.module,
            main=which.module, cex.main=2)
    par(mar=c(2,2.3,0.5,0.8))
    barplot(ME, col=which.module, main="", cex.main=2,
            ylab="eigengene expression",xlab="array sample")

}




########################################***Cell Chat Analysis***########################################
