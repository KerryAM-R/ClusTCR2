#' Code for plotting the Motif based on a specific CDR3 length and V gene (see [netplot] for )
#' @param ClusTCR Matrix file produce from [ClusTCR_list_to_matrix]
#' @param mcl_cluster matrix file produce from [mcl_to_lable_clusters]
#' @param Clust_column_name Name of clustering column from mcl_cluster file e.g. cluster
#' @param Clust_selected Select which cluster to label.
#' @param selected_col Color of selected cluster (Default = purple)
#' @param non_selected_col Color of selected cluster (Default = grey80)
#' @param selected_text_col Color of selected cluster text (Default = black)
#' @param non_selected_text_col Color of selected clusters text (Default = grey40)
#' @param selected_text_size Text size of selected cluster (Default = 3)
#' @param non_selected_text_size Text size of non-selected clusters (Default = 2)
#' @param label Name to display on cluster: Name (CDR3_V_gene_Cluster), cluster, CDR3, V_gene, Len (length of CDR3 sequence)
#' @param select_alpha Transparency of selected cluster
#' @param select_alpha Transparency of non-selected clusters
#' @importFrom  VLF aa.count.function
#' @importFrom motifStack pcm2pfm
#' @importFrom ggseqlogo ggseqlogo
#' @import ggplot2
#' @export

motif_plot <- function(ClusTCR, Clust_column_name="cluster",Clust_selected=NULL) {
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
