library(ggplot2)
library(GEOquery)
library(ggrepel)
library(ggpubr)
library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(reticulate)
library(MAESTRO)
library(harmony)
library(tidyverse)
library(cowplot)
library(CellChat)
library(NMF)
library(ggalluvial)
library(svglite)
library(ComplexHeatmap)
library(stringr)
c54 <- c('dodgerblue2','#FF7F00',
         '#FB9A99','#CAB2D6','khaki2','deeppink1','lightsalmon',      
         'steelblue4','yellow4','yellow3',
         'red2','orange','cornflowerblue', 'magenta','darkolivegreen4',
         'indianred1','tan4','darkblue','mediumorchid1','firebrick4',
         'tan3','tan1','darkgray',
         'wheat4','#DDAD4B','seagreen1','moccasin',
         'mediumvioletred','seagreen','cadetblue1','tan2',
         'tomato3','#7CE3D8','gainsboro','gold1','skyblue2',
         '#FDBF6F','gray70','darkorange4','orchid1',
         'darkturquoise','maroon','brown')
##Pre-process scRNA-seq data
GEO <- as.character(Dataset_information[1,"Dataset_ID"])
Datatype <- as.character(Dataset_information[1,"Data type"])
setwd(paste0("~/OralDB/GEO/ST/",GEO,"/"))
#Read multiple samples and create seurat objects
samples <- list.files('./sample_file/')
seurat_list <- list()
for (sample in samples){
  data.path <- paste0("./sample_file/",sample)  
  seurat_data <- Read10X(data.dir = data.path)
  seurat_obj <- CreateSeuratObject(counts=seurat_data,
                                   project = sample,
                                   min.features = 200,
                                   min.cells = 3)
  seurat_list <- append(seurat_list,seurat_obj)
}
seurat_metadata <- merge(seurat_list[[1]],
                         y = seurat_list[-1],
                         add.cell.ids = samples)
#QC
seurat_metadata[["percent.mt"]] <- PercentageFeatureSet(seurat_metadata,pattern = "^MT-")
VlnPlot(seurat_metadata,features = c("nCount_RNA","nFeature_RNA","percent.mt"),raster = TRUE)
seurat_metadata <- subset(seurat_metadata,percent.mt<25&nFeature_RNA>200&nFeature_RNA<10000)

#Extract the expression matrix
count <- GetAssayData(seurat_metadata,layer = "counts")
count@factors$metadata <- seurat_metadata@meta.data
save(count,file = "./count_GSE161267.Rdata")
#Standardization, PCA, batch correction, dimensionality reduction, clustering
seurat_metadata <- seurat_metadata %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures=4000) %>% 
  ScaleData(vars.to.regress="percent.mt") %>% 
  RunPCA(npcs=50)
ElbowPlot(seurat_metadata,ndims=50)
seurat_metadata <- seurat_metadata %>% 
  RunHarmony(dims=1:20,group.by.vars = "patient") %>% 
  RunUMAP(dims=1:20,reduction="harmony") %>% 
  FindNeighbors(reduction='harmony',dims=1:20) %>% 
  FindClusters(resolution=seq(0.1,2,0.1)) %>% 
  identity()
DimPlot(seurat_metadata,reduction="umap",label=T)+theme_bw()
#Integrate multiple data by JoinLayers()
seurat_metadata <- JoinLayers(seurat_metadata)
#Search for differential genes (manual annotation of cell types)
DEGs <- FindAllMarkers(seurat_metadata,logfc.threshold = 0.25,min.pct = 0.25,assay = "RNA")

Idents(seurat_metadata) <- "seurat_clusters"
seurat_metadata <- RenameIdents(seurat_metadata,"0"="Myofibroblasts",
                                "1"="Endothelial",
                                "2"="CD4_T",
                                "3"="Plasma",
                                "4"="CD8_T",
                                "5"="Myofibroblasts",
                                "6"="Mono_Macro",
                                "7"="Epithelial",
                                "8"="Endothelial",
                                "9"="Endothelial",
                                "10"="Endothelial",
                                "11"="B_cell",
                                "12"="Mast",
                                "14"="Myofibroblasts",
                                "15"="NK",
                                "16"="Epithelial",
                                "17"="Mono_Macro",
                                "18"="Epithelial",
                                "19"="Myofibroblasts",
                                "20"="Myofibroblasts",
                                "21"="Endothelial",
                                "22"="Plasma",
                                "23"="Mono_Macro",
                                "24"="Mono_Macro",
                                "26"="Myofibroblasts",
                                "27"="CD4_T",
                                "28"="Plasma"
)
seurat_metadata@meta.data$minor <- Idents(seurat_metadata)
save(seurat_metadata,file = paste0('./seurat_metadata_',GEO,'.Rdata'))

##Overview preprocess
#Order celltype 
order_name <- c("Epithelial", "Keratinocytes", "Salivary", "Fibroblasts", "Myofibroblasts", "MSCs","VSMCs", "Endothelial", "T_cells", "CD8+T_cells", "CD4+T_cells", "Treg", "Tprolif", "B_cells", "Plasma", "NK", "Myeloid", "Mono_Macro", "Monocytes", "Macrophages", "DCs", "pDCs", "mDCs", "Neutrophils", "Mast", "Myocytes", "Mural", "Astrocyte", "Schwann_cells", "LECs", "Progenitor", "Melanocytes", "Platelets", "others")
seurat_metadata$celltype <- seurat_metadata$minor
unique(seurat_metadata$celltype)
cell_types <- seurat_metadata@meta.data$celltype
common_types <- intersect(order_name, unique(cell_types))
cell_types <- factor(cell_types,levels = common_types)
seurat_metadata@meta.data$celltype <- cell_types
DimPlot(seurat_metadata,group.by = "celltype",label = T,raster = FALSE)
Idents(seurat_metadata) <- "celltype"
DimPlot(seurat_metadata,label = T,raster = FALSE)

#Order condition
table(seurat_metadata$condition)
seurat_metadata$condition <- factor(seurat_metadata$condition,levels = mixedsort(unique(seurat_metadata$condition))) 

#Order 
table(seurat_metadata$sample)
seurat_metadata$sample <- factor(seurat_metadata$sample,levels = mixedsort(unique(seurat_metadata$sample)))

#celltype data
celltype <- as.data.frame(row.names = rownames(seurat_metadata@meta.data),seurat_metadata@meta.data$celltype)
colnames(celltype) <- "celltype"
save(celltype,file = paste0("celltype_",GEO,".Rdata"))

#cluster_data for DEGs and Gene_corr
for(i in 1:nrow(Dataset_information)){
  GEO <- as.character(Dataset_information[i,"Dataset_ID"])
  Type <- as.character(Dataset_information[i,"Data type"])
  if(Type=="scRNA-seq"){
    setwd(paste0("~/OralDB/GEO/scRNA-seq/",GEO))
    load(paste0("./seurat_metadata_",GEO,".Rdata"))
    #get counts matrix
    count <- GetAssayData(seurat_metadata,layer = "counts")
    count@factors$metadata <- seurat_metadata@meta.data
    save(count,file = paste0("./count_",GEO,".Rdata"))
    #get metadata
    metadata <- seurat_metadata@meta.data
    save(metadata,file = paste0("./metadata_",GEO,".Rdata"))
    #get cluster data
    if(!file.exists("./cluster_data")){
      dir.create("./cluster_data")
    }
    load(paste0("./celltype_",GEO,".Rdata"))
    cluster <- as.character(unique(celltype$celltype))
    setwd("./cluster_data")
    for(c in cluster){
      cluster_celltype <- subset(seurat_metadata,idents = c)
      save(cluster_celltype,file = paste0("./cluster_celltype_",c,"_",GEO,".Rdata"))
    }
  }
}

#Overview Fig
load("~/OralDB/GEO/Dataset_information.Rdata")
for(i in 1:nrow(Dataset_information)){
  GEO <- as.character(Dataset_information[i,"Dataset name"])
  Type <- as.character(Dataset_information[i,"Data type"])
  if(Type == "scRNA-seq"){
    setwd(paste0("~/OralDB/GEO/",Type,"/",GEO))
    load(paste0("./seurat_metadata_",GEO,".Rdata"))
    celltype <- c("minor")
    for(ct in celltype){
      if(ct == "minor"){
        #minor_UMAP
        #order_name <- c("others","Mast","Endothelial","pDCs","Mono_Macro","NK","Plasma","B_cell","CD8_T","Treg","CD4_T","Epithelial","Myofibroblasts","Fibroblasts")
        p2 <- DimPlot(seurat_metadata,label = T,reduction = "umap",cols = c54,group.by = ct,raster = FALSE,label.size = 4.5,
                      order = order_name,repel = TRUE)+
          labs(title = paste0(GEO,"_",ct))+
          theme_bw()+
          xlab("UMAP1")+
          ylab("UMAP2")+
          theme(plot.title = element_text(size = 12,hjust = 0.5))
        #minor_sample_prop
        sample.prop <- as.data.frame(prop.table(table(seurat_metadata$minor,seurat_metadata$sample)))
        colnames(sample.prop) <- c("celltype","sample","proportion")
        p4 <- ggplot(sample.prop,aes(x=sample,y=proportion,fill=celltype))+
          geom_bar(stat = "identity",position = "fill")+
          ggtitle(paste0(GEO,"_minor_sample"))+
          theme_bw()+
          ylab("Proportion")+
          guides(fill=guide_legend(title = NULL))+
          scale_fill_manual(values = c54,limits = rev(intersect(order_name,unique(sample.prop$celltype))))+
          theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.ticks.length = unit(0.5,"cm"),plot.title = element_text(hjust = 0.5,size = 15))
        #minor_cell_number
        minor_cell_number <- as.data.frame(table(seurat_metadata$minor))
        colnames(minor_cell_number) <- c("celltype","number")
        minor_cell_number$celltype <- factor(minor_cell_number$celltype, levels = rev(intersect(order_name,unique(minor_cell_number$celltype))))
        minor_cell_number <- minor_cell_number[order(factor(minor_cell_number$celltype)),]
      }
    }
    ggsave(plot = p2,filename = paste0("./minor_UMAP_",GEO,".png"),width = 8,height = 6,dpi = 600)
    plot_grid(p1, p2, ncol = 2,align = "h")
    ggsave(filename = paste0("./","UMAP_",GEO,".png"),width = 14,height = 5,dpi = 600)
    plot_grid(p3,p4,ncol = 2,align = "h")
    ggsave(filename = paste0("./","sample_barplot_",GEO,".png"),width = 15,height = 5,dpi = 600)
    png(filename = paste0("./pieplot_",GEO,".png"),width = 3500,height = 1600,res = 200)
    par(mfrow=c(1,2),xpd=TRUE)
    pie(major_cell_number$number,labels = with(major_cell_number,paste0(celltype,"(",number,")")),col = c54,radius = 0.8,clockwise = T,init.angle = 1,main = "Major",cex = 1.2)+
      theme(
        plot.margin = margin(t = 0,  
                             r = 0,  
                             b = 0,  
                             l = 0,  
                             unit = "cm"))
    pie(minor_cell_number$number,labels = with(minor_cell_number,paste0(celltype,"(",number,")")),col = c54,radius = 0.8,clockwise = T,init.angle = 1,main = "Minor",cex = 1.2)+
      theme(
        plot.margin = margin(t = 0,
                             r = 0,
                             b = 0,
                             l = 0, 
                             unit = "cm"))
    dev.off()
  }
}

#Dotplot for marker gene
#Modify the dotplot function
DotPlot1 <- function (object, features, assay = NULL, cols = c("lightgrey","#9C65B2","#BC392D", 
                                                               "#511B66"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                      idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                      scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, cells = colnames(object[[assay]]), 
                                        idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- FetchData(object = object, vars = split.by)[cells, 
                                                          split.by]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop(paste0("Need to specify at least ", length(x = unique(x = splits)), 
                    " colors using the cols parameter"))
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = log1p(data.use))
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- unlist(x = lapply(X = data.plot$id, FUN = function(x) sub(paste0(".*_(", 
                                                                                   paste(sort(unique(x = splits), decreasing = TRUE), 
                                                                                         collapse = "|"), ")$"), "\\1", x)))
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  data.plot$id <- factor(data.plot$id, levels = rev(levels(data.plot$id)))
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by))  +theme_bw()+scale.func(range = c(0, dot.scale), 
                                                                                                                                                limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                                          axis.title.y = element_blank(),
                                                                                                                                                                                          panel.background = element_rect(fill = "white"),  plot.background = element_rect(fill= "white")   
                                                                                                                                                ) + guides(size = guide_legend(title = "Per. Exp(%)")) + 
    labs(x = "", y = ifelse(test = is.null(x = split.by), 
                            yes = "", no = "")) + theme_cowplot()+ theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))+theme(panel.grid.major = element_line(color = "#f9f9f9"), panel.grid.minor = element_line(color = "#f9f9f9"))
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                              space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
                                                                                           units = "lines"),axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Avg. Exp"))
  }
  return(plot)
}
Dataset_information <- subset(Dataset_information,Dataset_information$`Data type`=="scRNA-seq")
for(i in 1:nrow(Dataset_information)){
  marker <- c(
    "KRT14","KRT5","JUP","CD24","EPCAM","CDH1",
    "KRT7","SLPI","GSTA1","CLDN3",
    "CD24","KRT19",
    "DCN","COL1A1","CFD","PTGDS",
    "TAGLN","MYL9","ACTA2","S100A4","FN1","LUM",
    "NOTCH3","MCAM","FRZB",
    "PECAM1","VWF","CDH5",
    "CD3D","CD3E","CD8A","CD8B","GZMK","CCR7","LEF1","IL7R","CD4","CD27",
    "FOXP3","CTLA4","LAYN","IL2RA","IKZF2",
    "MKI67","TOP2A","CCND1",
    "MS4A1","CD19","CD79A",
    "MZB1","SDC1","XBP1","IGHG",
    "GNLY","NKG7","KLRC1",
    "LYZ","CD14","CD68","C1QC",
    "IL3RA","IRF4","SOX4","LILRA4",
    "HLA-DRA","HLA-DPB1","HLA-DQA1",
    "CXCL8","G0S2","CSF3R",
    "MS4A2","KIT","CPA3",
    "ACTA1","MYL1","MYH2",
    "CXCL14","CSRP2","SFRP1","MT2A",
    "MPZ","MBP","SOX10","NGR1",
    "PROX1","CCL21","LYVE1",
    "DCT","PMEL","TYRP1",
    "ITGA2B","ITGB3","GP1BA"
  )
  GEO <- as.character(Dataset_information[i,"Dataset_ID"])
  Type <- as.character(Dataset_information[i,"Data type"])
  setwd(paste0("~/OralDB/GEO/",Type,"/",GEO))
  load(paste0("./seurat_metadata_",GEO,".Rdata"))
  CairoPNG(paste0("~/OralDB/GEO/scRNA-seq/", GEO,"/Dotplot_",GEO,".png"), height = 1250, width = 3600, res = 300)
  DotPlot1(seurat_metadata,features = marker,scale = FALSE)+RotatedAxis()
  dev.off()
}

#DEGs_all
Dataset_information <- subset(Dataset_information,Dataset_information$`Data type`=="scRNA-seq")
for(i in 1:nrow(Dataset_information)){
  GEO <- as.character(Dataset_information[i,"Dataset_ID"])
  Type <- as.character(Dataset_information[i,"Data type"])
  path <- paste0("~/OralDB/GEO/",Type,"/",GEO)
  if(file.exists(path)){
    if(Type == "scRNA-seq"){
      setwd(paste0("~/OralDB/GEO/",Type,"/",GEO))
      load(paste0("./seurat_metadata_",GEO,".Rdata"))
      #Get gene expression matrix
      count <- GetAssayData(seurat_metadata,layer = "counts")
      count@factors$metadata <- seurat_metadata@meta.data
      save(count,file = paste0("./count_",GEO,".Rdata"))
      #DEGs
      DEGs <- FindAllMarkers(seurat_metadata,min.pct = 0.25,logfc.threshold = 0.25,assay = "RNA")
      DEGs <- DEGs[,c("gene","cluster","avg_log2FC","p_val","p_val_adj","pct.1","pct.2")]
      save(DEGs,file = paste0("./DEGs_",GEO,".Rdata"))
    }
  }
}

#DEGs_each_condition
for(i in 1:nrow(Dataset_information)){
  GEO <- as.character(Dataset_information1[i,"Dataset_ID"])
  Type <- as.character(Dataset_information1[i,"Data type"])
  if(Type=="scRNA-seq"){
    setwd(paste0("~/OralDB/GEO/scRNA-seq/",GEO))
    load(paste0("./metadata_",GEO,".Rdata"))
    condition = unique(metadata$condition)
    celltype = unique(metadata$celltype)
    if(!file.exists("./DEGs")){
      dir.create("./DEGs")
    }
    for(c in celltype){
      if(!file.exists(paste0("./DEGs/",c))){
        dir.create(paste0("./DEGs/",c))
      }
      load(paste0("./cluster_data/cluster_celltype_",c,"_",GEO,".Rdata"))
      for(i in 1:length(condition)){
        c1 = condition[i]
        for(j in 1:length(condition)){
          c2 = condition[j]
          if(c1 != c2){
            #Check the number of cells under each condition
            cells_c1 <- sum(cluster_celltype$condition == c1)
            cells_c2 <- sum(cluster_celltype$condition == c2)
            
            if(cells_c1 >= 3 & cells_c2 >= 3){
              tryCatch({
                DEGs <- FindMarkers(cluster_celltype, group.by = "condition", 
                                    ident.1 = c1, ident.2 = c2, 
                                    min.pct=0.25, logfc.threshold = 0.25)
                DEGs$gene <- rownames(DEGs)
                DEGs <- DEGs[which(abs(DEGs$avg_log2FC) >= 0.5 & DEGs$p_val_adj<0.05),]
                DEGs <- DEGs[,c('gene','p_val','p_val_adj','avg_log2FC','pct.1','pct.2')]
                DEGs$p_val <- format(DEGs$p_val, scientific = TRUE, digits = 2)
                DEGs$p_val_adj <- format(DEGs$p_val_adj, scientific = TRUE, digits = 2)
                DEGs$avg_log2FC <- round(DEGs$avg_log2FC,2)
                save(DEGs, file = paste0("./DEGs/",c,"/",c,"_",c1,"_vs_",c2,".Rdata"))
              }, error = function(e) {
                message(paste("Error in comparing", c1, "vs", c2, "for cell type", c, ":", e$message))
              })
            } else {
              message(paste("Skipping comparison: ",GEO, "_",c1, "vs", c2, "for cell type", c, 
                            "- insufficient cells (", cells_c1, "vs", cells_c2, ")"))
            }
          }
        }
      }
    }
  }
}

#GO Enrichment
for(i in 1:nrow(Dataset_information)){
  GEO <- as.character(Dataset_information[i,"Dataset_ID"])
  Type <- as.character(Dataset_information[i,"Data type"])
  if(Type=="scRNA-seq"){
    setwd(paste0("~/OralDB/GEO/scRNA-seq/",GEO))
    load(paste0("./metadata_",GEO,".Rdata"))
    condition = unique(metadata$condition)
    celltype = unique(metadata$celltype)
    if(!file.exists("./GO")){
      dir.create("./GO")
    }
    for(c in celltype){
      if(!file.exists(paste0("./GO/",c))){
        dir.create(paste0("./GO/",c))
      }
      for(i in 1:length(condition)){
        c1 = condition[i]
        for(j in 1:length(condition)){
          c2 = condition[j]
          if(c1 != c2){
            filename = paste0("./DEGs/",c,"/",c,"_",c1,"_vs_",c2,".Rdata")
            if(file.exists(filename)){
              load(filename)
              DEGs$p_val_adj <- as.numeric(DEGs$p_val_adj)
              
              DEGs_up <- DEGs[which(DEGs$avg_log2FC >= 1 & DEGs$p_val_adj<0.05),]
              DEGs_up <- DEGs_up[order(DEGs_up$avg_log2FC,decreasing = TRUE),]
              ego_bp <- data.frame()
              ego_BP <- enrichGO(gene          = row.names(DEGs_up),
                                 OrgDb         = "org.Hs.eg.db",
                                 keyType       = 'SYMBOL',
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05) 
              ego_bp <- data.frame(ego_BP)
              ego_bp_UP <- ego_bp
              DEGs_down <- DEGs[which(DEGs$avg_log2FC < -1 & DEGs$p_val_adj<0.05),]
              DEGs_down <- DEGs_down[order(DEGs_down$avg_log2FC,decreasing = TRUE),]
              ego_bp <- data.frame()
              ego_BP <- enrichGO(gene          = row.names(DEGs_down),
                                 OrgDb         = "org.Hs.eg.db",
                                 keyType       = 'SYMBOL',
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05) 
              ego_bp <- data.frame(ego_BP)
              ego_bp_DOWN <- ego_bp
            }else{
              ego_bp_UP <- data.frame()
              ego_bp_DOWN <- data.frame()
            }
            save(ego_bp_UP,file = paste0("./GO/",c,"/",c,"_",c1,"_vs_",c2,"_UP.Rdata"))
            save(ego_bp_DOWN,file = paste0("./GO/",c,"/",c,"_",c1,"_vs_",c2,"_DOWN.Rdata"))
          }
        }
      }
    }
  }
}

##Cell-cell interactions by Cellchat
#run CCI-----------
for(i in 1:nrow(Dataset_information)){
  GEO <- as.character(Dataset_information[i,"Dataset_ID"])
  Type <- as.character(Dataset_information[i,"Data type"])
  if(Type == "scRNA-seq"){
    setwd(paste0("~/OralDB/GEO/",Type,"/",GEO))
    load(paste0("./seurat_metadata_",GEO,".Rdata"))
    if (!file.exists("./CCI")) {
      dir.create("./CCI")
    }
    setwd("./CCI")
    pre_CCI()
    CCI_heatmap(GEO,cellchat)
    CCI_number(GEO,mat)
    CCI_bubbleplot(Type,GEO,cellchat)
  }
}
#pre-process-------------
pre_CCI <- function(){
  cellchat <- createCellChat(seurat_metadata,group.by = "celltype")
  groupSize <- as.numeric(table(cellchat@idents))
  cellchatDB <- CellChatDB.human
  unique(cellchatDB$interaction$annotation)
  cellchatDB.use <- subsetDB(cellchatDB,search = "Secreted Signaling")
  cellchat@DB <- cellchatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat,PPI.human)
  #Infer the cell communication network
  cellchat <- computeCommunProb(cellchat,raw.use = FALSE,population.size = TRUE)
  cellchat <- filterCommunication(cellchat,min.cells = 10)
  #df.net <- subsetCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  #df.netp <- subsetCommunication(cellchat,slot.name = "netP")
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  mat <- cellchat@net$count
}
##Interaction Count Heatmap----------------
CCI_heatmap <- function(GEO,cellchat){
  width = 110*12+1000
  height = 60*12+1200
  data <- as.data.frame(cellchat@net$count)
  png(filename = paste0("./Interaction_Heatmap_",GEO,".png"),width = width,height = height,res = 120)
  pdf("article/Interaction_Counts_Heatmap.pdf", width = 8, height = 7.5)
  Heatmap(cellchat@net[["count"]],name = "count",
          col = c("#2171B5","#F2F2F2", "#ED645A"),
          rect_gp =gpar(col="black",lwd =1),
          column_title_gp = gpar(fontsize = 15),
          column_title = "Interaction Counts",
          cell_fun = function(j, i, x, y, w, h, col){
            grid.text(data[i, j], x, y)
          }
  )
  dev.off()
}
##Single (number/weight of interaction)--------------
CCI_number <- function(GEO,mat){
  for(i in 1:nrow(mat)){
    mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
    mat2[i,] <- mat[i,]
    name <- rownames(mat)[i]
    if(grepl("\\+", name)){
      name <- gsub("\\+", "", name)
    }
    myfilename <- paste0(name,"_",GEO,".png")
    width = 2000
    height = 1800
    png(filename = myfilename,width = width,height = height,res = 300)
    netVisual_circle(mat2,
                     vertex.weight = groupSize,weight.scale = T,
                     arrow.width = 0.2,arrow.size = 0.1,
                     edge.weight.max = max(mat)
    )
    grid.text(rownames(mat)[i],x = unit(0.5, "npc"), y = unit(0.95, "npc"),gp = gpar(fontsize = 15))
    dev.off()
  }
}
##Target and Source Bubble Plot-----------------
CCI_bubbleplot <-function(Type,GEO,cellchat){
  if (!file.exists("./GSEA")) {
    dir.create("./bubble_plot")
  }
  setwd("./bubble_plot")
  type_list <- c("source","target")
  for(type in type_list){
    setwd(paste0("~/OralDB/GEO/",Type,"/",GEO,"/CCI"))
    src = paste0("~/OralDB/GEO/",Type,"/",GEO,"/CCI/bubble_plot/",type)
    if(!file.exists(src)){
      dir.create(src)
    }
    setwd(src)
    for(i in 1:length(levels(cellchat@idents))){
      use <- as.character()
      use <- levels(cellchat@idents)[-i]
      if(type == "source"){
        source_use <- levels(cellchat@idents)[i]
        target_use <- use
      }else{
        
        source_use <- use
        target_use <- levels(cellchat@idents)[i]
      }
      name <- levels(cellchat@idents)[i]
      if(grepl("\\+", name)){
        name <- gsub("\\+", "", name)
      }
      filename = paste0(type,"_",name,".png")
      return_data <- netVisual_bubble(cellchat,
                                      sources.use = source_use,
                                      targets.use = target_use,
                                      angle.x = 45,
                                      remove.isolate = FALSE,
                                      thresh = 0.01,
                                      title.name = str_to_title(type),
                                      font.size = 11,
                                      font.size.title = 15,
                                      return.data = TRUE
      )
      cci_row <- length(unique(return_data$communication$source.target))
      cci_col <- length(unique(return_data$communication$interaction_name_2))
      width <- cci_row*0.8
      height <- cci_col*0.25
      netVisual_bubble(cellchat,
                       sources.use = source_use,
                       targets.use = target_use,
                       angle.x = 45,
                       remove.isolate = FALSE,
                       thresh = 0.01,
                       title.name = str_to_title(type),
                       font.size = 11,
                       font.size.title = 15
      )
      ggsave(file = filename,width = width,height = height,dpi = 300,limitsize = FALSE)
    }
    
  }
} 


