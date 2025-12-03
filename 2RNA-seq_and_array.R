library(DESeq2)
library(stats)
library(affy)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)

##RNA-seq
#Count to TPM
CountData$gene <- rownames(CountData)
tmp <- inner_join(CountData,gene_length,by = 'gene')
#tmp <- count
TPM <- data.frame(row.names = tmp$gene)
for(i in 1:(dim(tmp)[2]-2)){
  col <- tmp[[i]]
  len <- tmp[[dim(tmp)[2]]]
  rate <- col/len
  N <- sum(col)
  TPM_i <- (rate*1e6)/(sum(rate))
  TPM_i <- pmax(TPM_i,0) %>% as.data.frame()
  colnames(TPM_i) <- colnames(tmp)[i]
  TPM <- cbind(TPM,TPM_i)
}
TPM <- as.matrix(TPM)
TPM <- Matrix::Matrix(TPM,sparse = T)
TPM@factors$condition <- colData
save(TPM,file = "./TPM_RNAseq.Rdata")

#DEGs
for(i in 1:nrow(Dataset_information)){
  GEO <- as.character(Dataset_information[i,"Dataset_ID"])
  Type <- as.character(Dataset_information[i,"Data type"])
  if(Type=="RNA-seq"){
    setwd(paste0("/home/cuihao/OralDB/GEO/RNA-seq/",GEO))
    load(paste0("./count_",GEO,".Rdata"))
    condition = unique(count@factors$condition$condition)
    condition <- as.character(condition)
    if(!file.exists("./DEGs")){
      dir.create("./DEGs")
    }
    load(paste0("./dds","_",GEO,".Rdata"))
    for(i in 1:length(condition)){
      c1 = condition[i]
      for(j in 1:length(condition)){
        c2 = condition[j]
        if(c1 != c2){
          res <- results(dds,contrast = c("condition",c1,c2))
          DEGs <- data.frame(res,stringsAsFactors = FALSE,check.names = FALSE)
          if(nrow(DEGs)>0){
            DEGs <- DEGs[order(DEGs$pvalue,DEGs$log2FoldChange,decreasing = c(FALSE,TRUE)),]
            DEGs$gene <- rownames(DEGs)
            DEGs <- DEGs[which(abs(DEGs$log2FoldChange) >= 0.5 & DEGs$padj<0.05),]
            DEGs <- DEGs[,c("gene","pvalue","padj","log2FoldChange","baseMean","stat","lfcSE")]
            #DEGs$pvalue <- format(DEGs$pvalue, scientific = TRUE, digits = 2)
            #DEGs$padj <- format(DEGs$padj, scientific = TRUE, digits = 2)
            DEGs$log2FoldChange <- round(DEGs$log2FoldChange,2)
            save(DEGs, file = paste0("./DEGs/",c1,"_vs_",c2,".Rdata"))
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
  if(Type=="RNA-seq"){
    setwd(paste0("/home/cuihao/OralDB/GEO/RNA-seq/",GEO))
    load(paste0("./count_",GEO,".Rdata"))
    condition = unique(count@factors$condition$condition)
    condition <- as.character(condition)
    if(!file.exists("./GO")){
      dir.create("./GO")
    }
    for(i in 1:length(condition)){
      c1 = condition[i]
      for(j in 1:length(condition)){
        c2 = condition[j]
        if(c1 != c2){
          filename = paste0("./DEGs/",c1,"_vs_",c2,".Rdata")
          if(file.exists(filename)){
            load(filename)
            DEGs_up <- DEGs[which(DEGs$log2FoldChange > 0 & DEGs$padj<0.05),]
            DEGs_up <- DEGs_up[order(DEGs_up$log2FoldChange,decreasing = TRUE),]
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
            DEGs_down <- DEGs[which(DEGs$log2FoldChange < 0 & DEGs$padj<0.05),]
            DEGs_down <- DEGs_down[order(DEGs_down$log2FoldChange,decreasing = TRUE),]
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
          save(ego_bp_UP,file = paste0("./GO/",c1,"_vs_",c2,"_UP.Rdata"))
          save(ego_bp_DOWN,file = paste0("./GO/",c1,"_vs_",c2,"_DOWN.Rdata"))
        }
      }
    }
  }
}

##Microarray
data.raw <- ReadAffy(celfile.path = "./raw_data/")
data.rma <- rma(data.raw)
expr <- exprs(data.rma)
expr <- as.data.frame(expr)
#Corresponding gene Symbol
expr$ID <- rownames(expr)
GPL570 <- getGEO("GPL570",destdir=".")
anno=GPL570@dataTable@table
probeSymbol<- data.frame(anno$ID,anno$`Gene Symbol`)
colnames(probeSymbol) <- c("ID","Symbol")
tmp=unlist(lapply(1:nrow(anno),function(i){
  tmp1=strsplit(anno[,"gene_assignment"][i],"///") #The splitting of multiple genes with one probe
  tmp2=unique(trimws(stringr::str_split(unlist(tmp1),"//",simplify = T)[,2])) #Extract the second symbol ID
  names(tmp2)=rep(anno[,1][i],length(tmp2))
  return(tmp2)
}))
probe2symbol=data.frame(ID=names(tmp),Symbol=tmp)
count <- merge(x = probe2symbol, y = expr,by = "ID")
count <- dplyr::select(count,-ID)
count <- as.data.frame(count)
count <- as.matrix(count)
count <- Matrix(count,sparse = T)
count@factors$condition <- colData
save(count,file = "count_array.Rdata")

#DEGs
for(i in 1:nrow(Dataset_information)){
  GEO <- as.character(Dataset_information[i,"Dataset_ID"])
  Type <- as.character(Dataset_information[i,"Data type"])
  if(Type=="array"){
    setwd(paste0("/home/cuihao/OralDB/GEO/",Type,"/",GEO))
    load(paste0("./count_",GEO,".Rdata"))
    condition = unique(count@factors$condition$condition)
    condition <- as.character(condition)
    if(!file.exists("./DEGs")){
      dir.create("./DEGs")
    }
    load(paste0("./dds","_",GEO,".Rdata"))
    for(i in 1:length(condition)){
      c1 = condition[i]
      for(j in 1:length(condition)){
        c2 = condition[j]
        if(c1 != c2){
          res <- results(dds,contrast = c("condition",c1,c2))
          DEGs <- data.frame(res,stringsAsFactors = FALSE,check.names = FALSE)
          if(nrow(DEGs)>0){
            DEGs <- DEGs[order(DEGs$pvalue,DEGs$log2FoldChange,decreasing = c(FALSE,TRUE)),]
            DEGs$gene <- rownames(DEGs)
            DEGs <- DEGs[which(abs(DEGs$log2FoldChange) >= 0.5 & DEGs$padj<0.05),]
            DEGs <- DEGs[,c("gene","pvalue","padj","log2FoldChange","baseMean","stat","lfcSE")]
            #DEGs$pvalue <- format(DEGs$pvalue, scientific = TRUE, digits = 2)
            #DEGs$padj <- format(DEGs$padj, scientific = TRUE, digits = 2)
            DEGs$log2FoldChange <- round(DEGs$log2FoldChange,2)
            save(DEGs, file = paste0("./DEGs/",c1,"_vs_",c2,".Rdata"))
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
  if(Type=="array"){
    setwd(paste0("/home/cuihao/OralDB/GEO/",Type,"/",GEO))
    load(paste0("./count_",GEO,".Rdata"))
    condition = unique(count@factors$condition$condition)
    condition <- as.character(condition)
    if(!file.exists("./GO")){
      dir.create("./GO")
    }
    for(i in 1:length(condition)){
      c1 = condition[i]
      for(j in 1:length(condition)){
        c2 = condition[j]
        if(c1 != c2){
          filename = paste0("./DEGs/",c1,"_vs_",c2,".Rdata")
          if(file.exists(filename)){
            load(filename)
            DEGs_up <- DEGs[which(DEGs$log2FoldChange > 0 & DEGs$padj<0.05),]
            DEGs_up <- DEGs_up[order(DEGs_up$log2FoldChange,decreasing = TRUE),]
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
            DEGs_down <- DEGs[which(DEGs$log2FoldChange < 0 & DEGs$padj<0.05),]
            DEGs_down <- DEGs_down[order(DEGs_down$log2FoldChange,decreasing = TRUE),]
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
          save(ego_bp_UP,file = paste0("./GO/",c1,"_vs_",c2,"_UP.Rdata"))
          save(ego_bp_DOWN,file = paste0("./GO/",c1,"_vs_",c2,"_DOWN.Rdata"))
        }
      }
    }
  }
}

##PCA analysis
Dataset_information <- subset(Dataset_information,Dataset_information$`Data type`!= "scRNA-seq")
for(i in 1:nrow(Dataset_information)){
  datatype <- Dataset_information[i,4]$`Data type`
  GEO <- Dataset_information[i,2]$`Dataset name`
  src = paste0("~/GEO/",datatype,"/",GEO)
  setwd(dir = src)
  if(GEO != ""){
    file_dds <- paste0("./dds_",GEO,".Rdata")
    load(file = file_dds)
    file_count <- paste0("./count_",GEO,".Rdata")
    load(file = file_count)
    rld <- DESeq2::rlog(dds,blind = TRUE)
    rld_mat <- assay(rld)
    pca <- prcomp(t(rld_mat))
    pca_data <- as.data.frame(pca$x)
    condition <- count@factors$condition
    condition = factor(condition$condition)
    title <- paste0("PCA_",GEO)
    percentage <- round(pca$sdev/sum(pca$sdev)*100,2)
    percentage <- paste0(colnames(pca_data),"(",as.character(percentage),"%",")")
    ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
      geom_point() +
      theme_bw()+
      xlab(percentage[1])+
      ylab(percentage[2])+
      ggtitle(title)+
      stat_ellipse(level = 0.95,na.rm = TRUE)+
      theme(
        panel.background = element_rect(fill = "transparent"),  
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid = element_blank(),  
        legend.background = element_rect(fill = "transparent")
      )
    ggsave(filename = paste0("/PCA_",GEO,".png"),width = 6,height = 4.5,dpi = 300)
  }
}

##Gene_corr
GEO <- as.character(Dataset_information[i,"Dataset_ID"])
Datatype <- as.character(Dataset_information[i,"Data type"])
type = "TPM_"
i=1
filename <- paste0("GEO/",Datatype,"/",GEO,"/",type,GEO,".Rdata")
load(file = filename)
data <- data.frame(TPM[input$Gene1,],TPM[input$Gene2,])  #input$Gene1 and input$Gene2 from the web page UI selection input box
colnames(data) <- c(paste0(input$Gene1),paste0(input$Gene2))
data[input$Gene1] <- log2(data[input$Gene1]+1)
data[input$Gene2] <- log2(data[input$Gene2]+1)
p <- ggscatter(data, x = input$Gene1, y = input$Gene2,
               color = "black",
               conf.int = TRUE, 
               add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
               cor.coef = T,
               cor.coeff.args = list(),
               cor.method = "pearson",size = 2,cor.coef.size = 5,add = "reg.line"
)+theme_bw()+
  xlab(paste0("Expression of ",input$Gene1))+
  ylab(paste0("Expression of ",input$Gene2))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
