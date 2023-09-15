#' Code for displaying the network.
#' @param ClusTCR File produced from mcl_cluster
#' @param Clust_column_name Name of clustering column from mcl_cluster file e.g. cluster (Re-numbering the original_cluster), Original_cluster, Clust_size_order (Based on cluster size e.g. number of nodes)
#' @param Clust_selected Select which cluster to label.
#' @param selected_col Color of selected cluster (Default = purple)
#' @param non_selected_col Color of selected cluster (Default = grey80)
#' @param selected_text_col Color of selected cluster text (Default = black)
#' @param non_selected_text_col Color of selected clusters text (Default = grey40)
#' @param selected_text_size Text size of selected cluster (Default = 3)
#' @param non_selected_text_size Text size of non-selected clusters (Default = 2)
#' @param label Name to display on cluster: Name (CDR3_V_gene_Cluster), cluster, CDR3, V_gene, Len (length of CDR3 sequence),CDR3_selected,Name_selected,cluster_selected, (_selected only prints names of the chosen cluster), None
#' @param alpha_selected Transparency of selected cluster (default = 1)
#' @param alpha_non_selected Transparency of non-selected clusters (default = 0.5)
#' @param colour Colour selected = "color_test" or all = "color_all"
#' @param all.colour Colours all points by: rainbow, random, heat.colors, terrain.colors, topo.colors, hcl.colors and default
#' @param filter_plot Filter's plot to remove connects grater than # e.g. 2 = 3 or more connections.
#' @importFrom ggnet ggnet2
#' @import ggplot2
#' @importFrom stringr str_sub
#' @importFrom network as.network
#' @export

netplot_ClusTCR2 <- function(ClusTCR, filter_plot = 0, Clust_selected=1,selected_col="purple",selected_text_col="black",selected_text_size=3,non_selected_text_size=2, Clust_column_name="cluster", label = c("Name","cluster","CDR3","V_gene","Len"), non_selected_col="grey80",non_selected_text_col="grey40",alpha_selected=1,alpha_non_selected=0.5, colour = "color_test",all.colour = "default") {

  gg_fill_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  } # colouring function

  net.initial <- ClusTCR[[2]]
  df_clust.initial <- ClusTCR[[1]]
  names(ClusTCR)

  df_clust <- subset(df_clust.initial,df_clust.initial$count>filter_plot)
  net.initial2 <- (net.initial[,colnames(net.initial) %in% df_clust$CDR3_Vgene])
  net2 <- (net.initial2[rownames(net.initial2) %in% df_clust$CDR3_Vgene,])
  net2

  col_unique <- as.data.frame(unique(df_clust[,names(df_clust) %in% Clust_column_name]))
  col_unique <- as.data.frame(col_unique[order(col_unique[,1]),])
  names(col_unique) <- Clust_column_name

  if (all.colour == "rainbow") {
    col_unique$col <- rev(rainbow(dim(col_unique)[1]))
  }

 else if (all.colour == "heat.colors") { #(c("rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors"))
   col_unique$col <- heat.colors(dim(col_unique)[1])
 }

  else if (all.colour == "terrain.colors") { #(c("rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors"))
    col_unique$col <- col.terrain(dim(col_unique)[1])
  }

  else if (all.colour == "topo.colors") { #(c("rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors"))
    col_unique$col <- topo.colors(dim(col_unique)[1])
  }

  else if (all.colour == "hcl.colors") { #(c("rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors"))
    col_unique$col <- hcl.colors(dim(col_unique)[1], palette = "viridis")
  }

  else {
    col_unique$col <- gg_fill_hue(dim(col_unique)[1])
  }

  df_clust <- merge(df_clust,col_unique,by = Clust_column_name)
  df_clust <- df_clust[order(df_clust$order),]
  colnames(net2) <- paste(colnames(net2), df_clust[,names(df_clust) %in% c(Clust_column_name)], df_clust$col,sep = "_")



  net = as.network(net2,ignore.eval = FALSE, loops = F, names.eval = 'testValue')


  z <- net$gal$n
  for (i in 1:z) {
    net$val[[i]]$Name <- gsub(".* ", "", net$val[[i]]$vertex.names)
  }
  for (i in 1:z) {
    net$val[[i]]$Name <- str_sub(net$val[[i]]$Name,0,-9)
  }

  for (i in 1:z) {
    net$val[[i]]$CDR3 <- str_split_fixed(net$val[[i]]$vertex.names, '_', 3)[1]
  }

  for (i in 1:z) {
    net$val[[i]]$None <- ""
  }

  for (i in 1:z) {
    net$val[[i]]$Len <- nchar(str_split_fixed(net$val[[i]]$vertex.names, '_', 4)[1])
  }
  for (i in 1:z) {
    net$val[[i]]$cluster <- str_split_fixed(net$val[[i]]$vertex.names, '_', 4)[3]
  }
  for (i in 1:z) {
    net$val[[i]]$V_gene <- str_split_fixed(net$val[[i]]$vertex.names, '_', 4)[2]
  }

  for (i in 1:z) {

    net$val[[i]]$color_all <- str_split_fixed(net$val[[i]]$vertex.names, '_', 4)[4]
  }

  for (i in 1:z) {
    net$val[[i]]$color_test <- ifelse(net$val[[i]]$cluster %in% Clust_selected, selected_col, non_selected_col)
  }

  for (i in 1:z) {
    net$val[[i]]$Name_selected <- ifelse(net$val[[i]]$cluster %in% Clust_selected, net$val[[i]]$Name, "")
  }

  for (i in 1:z) {
    net$val[[i]]$CDR3_selected <- ifelse(net$val[[i]]$cluster %in% Clust_selected, net$val[[i]]$V_gene, "")
  }

  for (i in 1:z) {
    net$val[[i]]$cluster_selected <- ifelse(net$val[[i]]$cluster %in% Clust_selected, net$val[[i]]$cluster, "")
  }

  for (i in 1:z) {
    net$val[[i]]$label_test <- ifelse(net$val[[i]]$cluster %in% Clust_selected, selected_text_col, non_selected_text_col)
  }
  for (i in 1:z) {
    net$val[[i]]$alpha_node <- ifelse(net$val[[i]]$cluster %in% Clust_selected, alpha_selected, alpha_non_selected)
  }
  for (i in 1:z) {
    net$val[[i]]$lab_size <- ifelse(net$val[[i]]$cluster %in% Clust_selected, selected_text_size, non_selected_text_size)
  }
  for (i in 1:z) {
    net$val[[i]]$lab_selected <- ifelse(net$val[[i]]$cluster %in% Clust_selected, net$val[[i]]$Name, net$val[[i]]$cluster)
  }

  ggnet2(net, size = "degree", label = label[1], color = colour, label.color = "label_test", node.alpha = "alpha_node", label.size = "lab_size")

}
