#' Code for plotting the Motif based on a specific CDR3 length and V gene (see [netplot_ClusTCR2] for ).
#' @param ClusTCR Matrix file produce from [mcl_cluster]
#' @param Clust_column_name Name of clustering column from mcl_cluster file e.g. cluster
#' @param Clust_selected Select which cluster to display. Only one at a time.
#' @importFrom  VLF aa.count.function
#' @importFrom ggseqlogo ggseqlogo
#' @import ggplot2
#' @export

motif_plot <- function(ClusTCR, Clust_column_name="cluster",Clust_selected=NULL) {
  source(system.file("Functions","motifStack.functions.R",package = "ClusTCR2"))
  library(VLF)
  library(ggseqlogo)


  net2 <- ClusTCR[[2]]
  df_clust <- ClusTCR[[1]]
  colnames(net2) <- paste(colnames(net2),df_clust[,names(df_clust) %in% Clust_column_name],sep="_")

  motif_DF <- as.data.frame(unique(colnames(net2)))
  names(motif_DF) <- "V1"

  z1 <- as.data.frame(t(as.data.frame(strsplit(motif_DF$V1,"_"))))
  head(z1)
  names(z1) <- c("V1","V2","V3")

  motif_DF$motif <- z1$V1
  motif_DF$cluster <- z1$V3

  motif_DF$len1 <- nchar(motif_DF[,names(motif_DF) %in% "motif"])
  if (length(which(motif_DF$cluster==c(Clust_selected)))>0) {


  motif_DF_interest <- motif_DF[motif_DF$cluster %in% c(Clust_selected),]

  # subset(motif_DF,motif_DF$cluster==Clust_selected)
  motif <- as.data.frame(t(as.data.frame(strsplit(motif_DF_interest$motif, ""))))
  len <- nchar(unique(motif_DF_interest$motif))[1]

  motif_count <- aa.count.function(cbind(x=1,y=2,motif),len)

  motif_count1_aa<-pcm2pfm(motif_count)


  ggseqlogo(motif_count, seq_type='aa', method='p') +
    ylab('bits')+
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    theme(
      axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
      axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
      axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
      axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
      legend.title  =element_blank(),
      legend.position = "right",
      legend.text = element_text(colour="black", size=12,family="serif")
    )

  }

  else {
    paste0("Cluster, ", Clust_selected,", does not exist")

  }
}

#' Code for plotting the Motif based on a specific CDR3 length and V gene (see [netplot_ClusTCR2] for ).
#' @param ClusTCRFile_large Matrix file produce from [mcl_cluster_large]
#' @param Clust_column_name Name of clustering column from mcl_cluster file e.g. cluster
#' @param Clust_selected Select which cluster to display. Only one at a time.
#' @importFrom  VLF aa.count.function
#' @importFrom ggseqlogo ggseqlogo
#' @import ggplot2
#' @export

motif_plot_large <- function(ClusTCRFile_large, Clust_column_name="Clust_size_order",Clust_selected=NULL) {
  source(system.file("Functions","motifStack.functions.R",package = "ClusTCR2"))
  library(VLF)
  library(ggseqlogo)

  df_clust <- ClusTCRFile_large[[1]]
  df_clust$Select_cluster_Name <- df_clust[,names(df_clust) %in% Clust_column_name]

  motif_df <- subset(df_clust,df_clust$Select_cluster_Name == Clust_selected)

  motif_DF <- as.data.frame(unique(motif_df$CDR3_Vgene))
  names(motif_DF) <- "V1"

  z1 <- as.data.frame(t(as.data.frame(strsplit(motif_DF$V1,"_"))))
  names(z1) <- c("V1","V2")
  motif_DF$motif <- z1$V1
  motif_DF$V_gene <- z1$V2
  motif_DF$selected <- Clust_selected
  head(motif_DF)
  motif_DF$len1 <- nchar(motif_DF[,names(motif_DF) %in% "motif"])
  motif <- as.data.frame(t(as.data.frame(strsplit(motif_DF$motif, ""))))

  len <- nchar(unique(motif_DF$motif))[1]

  motif_count <- aa.count.function(cbind(x=1,y=2,motif),len)

  motif_count1_aa<-pcm2pfm(motif_count)

  title.name <- unique(motif_DF$V_gene )

  ggseqlogo(motif_count, seq_type='aa', method='p') +
    ylab('bits')+
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    theme(
      axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
      axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
      axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
      axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
      legend.title  =element_blank(),
      legend.position = "right",
      title = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
      legend.text = element_text(colour="black", size=12,family="serif")
    ) +
    ggtitle(title.name)
}
