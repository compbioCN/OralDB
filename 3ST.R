library(Seurat)
library(hdf5r)
library(ggplot2)
library(data.table)
library(dplyr)
library(sctransform)
library(tidydr)
library(CARD)
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
GEO <- as.character(Dataset_information[1,"Dataset_ID"])
Datatype <- as.character(Dataset_information[1,"Data type"])
setwd(paste0("~/OralDB/GEO/ST/",GEO,"/"))
sample_list <- list()
n = length(sample_list)
for(i in 1:n){
  image = Read10X_Image(paste0("./raw_data/patient",i,"/spatial/"),image.name = "tissue_hires_image.png",filter.matrix = TRUE)
  sample = Read10X(paste0("./raw_data/patient",i,"/filtered_feature_bc_matrix/"))
  sample <- CreateSeuratObject(counts = sample, assay = "Spatial")
  sample$orig.ident <- paste0("patient",i)
  sample_list <- append(sample_list,sample)
}
for(i in 1:n){
  sample_list[[i]] <- sample_list[[i]][,sample_list[[i]]$nCount_Spatial>=5&sample_list[[i]]$nCount_Spatial<=15000] 
  sample_list[[i]] <- SCTransform(sample_list[[i]],assay = "Spatial")
  sample_list[[i]]@meta.data$sample <- names(sample_list)[i]
}

sample_ST <- merge(sample_list[[1]],sample_list[[2]])
for(i in 3:n){
  sample_ST <- merge(sample_ST,sample_list[[i]])
}

ST <- PrepSCTFindMarkers(sample_ST,assay = "SCT")
ST <- RunPCA(ST,npcs = 50)
ElbowPlot(ST)
ST <- ST %>% 
  RunHarmony(group.by.vars = "orig.ident",assay.use = "SCT") %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution=seq(0.5,1.5,0.1)) %>% 
  RunUMAP(dims=1:15,verbose=T)
Idents(ST) <- "SCT_snn_res.1"
DimPlot(ST,label = T,cols = c54,pt.size = 0.1)+
  theme_dr()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
SpatialDimPlot(ST,pt.size.factor = 2.5,crop = TRUE,image.scale = "hires",images = "A1")
for(i in 1:n){
  sample_name <- unique(ST$sample)[i]
  ST_patient <- subset(ST,sample == sample_name)
  save(ST_patient,file = paste0("./ST_sample/",sample_name,".Rdata"))
}

#CARD deconvolution
for(i in 1:8){
  #st
  ST_sample <- subset(ST,orig.ident == names(ST@images[i]))
  ST_count <- ST_sample@assays$Spatial$counts
  colnames(ST_count) <- sub("_1.*","",colnames(ST_count))
  ST_loca <- ST_sample@images[[1]]$centroids@coords
  ST_loca <- data.frame(ST_loca)
  # ST_loca1 <- ST_loca
  # ST_loca$x <- ST_loca1$y
  # ST_loca$y <- -ST_loca1$x
  rownames(ST_loca) <- ST_sample@images[[1]]$centroids@cells
  rownames(ST_loca) <- sub("_1.*","",rownames(ST_loca))
  spots_loca <- rownames(ST_loca)
  spots_use  <- spots_loca[spots_loca %in% colnames(ST_count)]
  
  ST_count_aligned <- ST_count[, spots_use, drop = FALSE]
  ST_count_aligned <- ST_count_aligned[, match(spots_use, colnames(ST_count_aligned)), drop = FALSE] 
  ST_loca_aligned  <- ST_loca[spots_use, , drop = FALSE]
  ST_loca <- ST_loca_aligned
  ST_count <- ST_count_aligned
  print(paste0(names(ST@images[i]),":", length(ST_loca$x)))
  #sc
  sc_count <- seurat_metadata@assays$RNA$counts
  sc_meta <- seurat_metadata@meta.data %>%
    dplyr::select(orig.ident,celltype)
  #rownames(sc_meta) <- sub(".*_","",rownames(sc_meta))
  sc_meta$cell_id <- rownames(sc_meta)
  #CARD
  CARD_obj <- createCARDObject(sc_count = sc_count,
                               sc_meta = sc_meta,
                               spatial_count = ST_count,
                               spatial_location = ST_loca,
                               ct.varname = "celltype",
                               ct.select = unique(sc_meta$celltype),
                               sample.varname = "orig.ident")
  CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)
  filename <- paste0("./ST_deconvolution/deconvolution_data/CARD_",names(ST@images[i]),".Rdata")
  save(CARD_obj,file = filename)
}
for(i in 1:8){
  p <- SpatialDimPlot(ST,pt.size.factor = 3,crop = TRUE,image.scale = "hires",images = names(ST@images[i]))
  p <- suppressMessages(p+theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                panel.background = element_blank(),
                                plot.background = element_blank(),
                                legend.background = element_blank(),
                                axis.text =element_blank(),
                                axis.ticks =element_blank(),
                                axis.title =element_blank(),
                                legend.title=element_text(size = 10,face="bold",family = "Helvetica"),
                                legend.text=element_text(size = 9,family = "Helvetica"),
                                strip.text = element_text(size = 16,face="bold"),
                                legend.position="right"))
  ggsave(p,filename = paste0("./ST_cluster/",names(ST@images[i]),"_cluster.png"),width = 6,height = 5,dpi = 600)
}

for(i in 1:8){
  filename = paste0("./ST_deconvolution/deconvolution_data/CARD_",names(ST@images[i]),".Rdata")
  load(filename)
  p1 <- CARD.visualize.pie2(proportion = CARD_obj@Proportion_CARD,spatial_location = CARD_obj@spatial_location,colors = c54,round_size = 18)
  p1
  ggsave(p1,filename = paste0("./ST_deconvolution/deconvolution_plot/",names(ST@images[i]),"_deconvo.png"),width = 8,height = 6,dpi = 600)
}
library(gtools)
library(scatterpie)
#Modify the CARD.visualize.pie2() function
CARD.visualize.pie2 <- function(proportion,spatial_location,colors = NULL,round_size){
  res_CARD = as.data.frame(proportion)
  res_CARD = res_CARD[,mixedsort(colnames(res_CARD))]
  location = as.data.frame(spatial_location)
  if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorCandidate = c("#1e77b4","#ff7d0b","#ceaaa3","#2c9f2c","#babc22","#d52828","#9267bc",
                     "#8b544c","#e277c1","#d42728","#adc6e8","#97df89","#fe9795","#4381bd","#f2941f","#5aa43a","#cc4d2e","#9f83c8","#91675a",
                     "#da8ec8","#929292","#c3c237","#b4e0ea","#bacceb","#f7c685",
                     "#dcf0d0","#f4a99f","#c8bad8",
                     "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
                     "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
                     "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
                     "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785")
  if(is.null(colors)){
    #colors = brewer.pal(11, "Spectral")
    if(ncol(res_CARD) > length(colorCandidate)){
      colors = colorRampPalette(colorCandidate)(ncol(res_CARD))
    }else{
      colors = colorCandidate[sample(1:length(colorCandidate),ncol(res_CARD))]
    }
  }else{
    colors = colors
  }
  data = cbind(res_CARD,location)
  data$x <- -data$x
  ct.select = colnames(res_CARD)
  p = suppressMessages(ggplot() + geom_scatterpie(aes(x=y, y=x,r = round_size),data=data, ###r默认是0.52
                                                  cols=ct.select,color=NA) + coord_fixed(ratio = 1) + 
                         scale_fill_manual(values =  colors)+
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               legend.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 16,face="bold"),
                               legend.text=element_text(size = 15),
                               legend.key = element_rect(colour = "grey89", fill = NA),
                               legend.key.size = unit(0.45, 'cm'),
                               strip.text = element_text(size = 16,face="bold"),
                               legend.position="right")+
                         guides(fill=guide_legend(title="Cell Type")))
  return(p)
}


