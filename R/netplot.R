#' Code for displaying the network.
#'
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
#' @import ggplot2
#' @import scales
#' @importFrom stringr str_sub
#' @import network
#' @import RColorBrewer
#' @import grDevices
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
    col_unique$col <- terrain.colors(dim(col_unique)[1])
  }

  else if (all.colour == "topo.colors") { #(c("rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors"))
    col_unique$col <- topo.colors(dim(col_unique)[1])
  }

  else if (all.colour == "hcl.colors") { #(c("rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors"))
    col_unique$col <- hcl.colors(dim(col_unique)[1], palette = "viridis")
  }
  else if (all.colour == "random") { #(c("rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors"))
    col_unique$col <- distinctColorPalette(dim(col_unique)[1])
  }
  else {
    col_unique$col <- gg_fill_hue(dim(col_unique)[1])
  }

  df_clust <- merge(df_clust,col_unique,by = Clust_column_name)
  df_clust <- df_clust[order(df_clust$order),]
  colnames(net2) <- paste(colnames(net2), df_clust[,names(df_clust) %in% c(Clust_column_name)], df_clust$col,sep = "_")

  net = network::as.network(net2,ignore.eval = FALSE, loops = F, names.eval = 'testValue')

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

#' Copied code from ggnet's ggnet2 function
#' @name ggnet2
#' @param net net plot from step 2.
#' @param mode Name of clustering column from mcl_cluster file e.g. cluster (Re-numbering the original_cluster), Original_cluster, Clust_size_order (Based on cluster size e.g. number of nodes)
#' @param mode = "fruchtermanreingold"
#' @param layout.par = NULL,
#' @param layout.exp = 0
#' @param layout.exp = 0
#' @param alpha = 1
#' @param color = "grey75"
#' @param shape = 19
#' @param size = 9
#' @param max_size = 9
#' @param na.rm = NA
#' @param palette = NULL
#' @param alpha.palette = NULL
#' @param alpha.legend = NA
#' @param color.palette = palette
#' @param color.legend = NA
#' @param shape.palette = NULL
#' @param shape.legend = NA
#' @param size.palette = NULL
#' @param size.legend = NA
#' @param size.zero = FALSE
#' @param size.cut = FALSE
#' @param size.min = NA
#' @param size.max = NA
#' @param label = FALSE
#' @param label.alpha = 1
#' @param label.color = "black"
#' @param label.size = max_size/2
#' @param label.trim = FALSE
#' @param node.alpha = alpha
#' @param node.color = color
#' @param node.label = label
#' @param node.shape = shape
#' @param node.size = size
#' @param edge.alpha = 1
#' @param edge.color = "grey50"
#' @param edge.lty = "solid"
#' @param edge.size = 0.25
#' @param edge.label = NULL
#' @param edge.label.alpha = 1
#' @param edge.label.color = label.color
#' @param edge.label.fill = "white"
#' @param edge.label.size = max_size/2
#' @param arrow.size = 0
#' @param arrow.gap = 0
#' @param arrow.type = "closed"
#' @param legend.size = 9
#' @param legend.position = "right"
#' @param ... Other functions in ggplot2
#' @import network
#' @importFrom sna degree gplot.layout.fruchtermanreingold
#' @import scales
#' @import RColorBrewer
#' @export

ggnet2 <- function (net, mode = "fruchtermanreingold", layout.par = NULL,
          layout.exp = 0, alpha = 1, color = "grey75", shape = 19,
          size = 9, max_size = 9, na.rm = NA, palette = NULL, alpha.palette = NULL,
          alpha.legend = NA, color.palette = palette, color.legend = NA,
          shape.palette = NULL, shape.legend = NA, size.palette = NULL,
          size.legend = NA, size.zero = FALSE, size.cut = FALSE, size.min = NA,
          size.max = NA, label = FALSE, label.alpha = 1, label.color = "black",
          label.size = max_size/2, label.trim = FALSE, node.alpha = alpha,
          node.color = color, node.label = label, node.shape = shape,
          node.size = size, edge.alpha = 1, edge.color = "grey50",
          edge.lty = "solid", edge.size = 0.25, edge.label = NULL,
          edge.label.alpha = 1, edge.label.color = label.color, edge.label.fill = "white",
          edge.label.size = max_size/2, arrow.size = 0, arrow.gap = 0,
          arrow.type = "closed", legend.size = 9, legend.position = "right",
          ...)
{

  net = try(network::network(net), silent = TRUE)

  get_v  <- function (x, attrname) {
    network::get.vertex.attribute(x, attrname = attrname)
  }

  get_e <- function (x, attrname) {
    network::get.edge.value(x, attrname = attrname)
  }
  set_mode = function(x, mode = network::get.network.attribute(x,
                                                               "bipartite")) {
    c(rep("actor", mode), rep("event", n_nodes - mode))
  }
  set_node = function(x, value, mode = TRUE) {
    if (is.null(x) || is.na(x) || is.infinite(x) || is.nan(x)) {
      stop(paste("incorrect", value, "value"))
    }
    else if (is.numeric(x) && any(x < 0)) {
      stop(paste("incorrect", value, "value"))
    }
    else if (length(x) == n_nodes) {
      x
    }
    else if (length(x) > 1) {
      stop(paste("incorrect", value, "length"))
    }
    else if (x %in% v_attr) {
      get_v(net, x)
    }
    else if (mode && x == "mode" & is_bip) {
      set_mode(net)
    }
    else {
      x
    }
  }
  set_edge = function(x, value) {
    if (is.null(x) || is.na(x) || is.infinite(x) || is.nan(x)) {
      stop(paste("incorrect", value, "value"))
    }
    else if (is.numeric(x) && any(x < 0)) {
      stop(paste("incorrect", value, "value"))
    }
    else if (length(x) == n_edges) {
      x
    }
    else if (length(x) > 1) {
      stop(paste("incorrect", value, "length"))
    }
    else if (x %in% e_attr) {
      get_e(net, x)
    }
    else {
      x
    }
  }
  set_attr = function(x) {
    if (length(x) == n_nodes) {
      x
    }
    else if (length(x) > 1) {
      stop(paste("incorrect coordinates length"))
    }
    else if (!x %in% v_attr) {
      stop(paste("vertex attribute", x, "was not found"))
    }
    else if (!is.numeric(get_v(net, x))) {
      stop(paste("vertex attribute", x, "is not numeric"))
    }
    else {
      get_v(net, x)
    }
  }
  set_name = function(x, y) {
    z = length(x) == 1 && x %in% v_attr
    z = ifelse(is.na(y), z, y)
    z = ifelse(isTRUE(z), x, z)
    ifelse(is.logical(z), "", z)
  }
  set_size = function(x) {
    y = x + (0 %in% x) * !size.zero
    y = scales::rescale_max(y)
    y = (scales::abs_area(max_size))(y)
    if (is.null(names(x)))
      names(y) = x
    else names(y) = names(x)
    y
  }
  is_one = function(x) length(unique(x)) == 1
  is_col = function(x) all(is.numeric(x)) | all(network::is.color(x))
  n_nodes = network::network.size(net)
  n_edges = network::network.edgecount(net)
  v_attr = network::list.vertex.attributes(net)
  e_attr = network::list.edge.attributes(net)
  is_bip = network::is.bipartite(net)
  is_dir = ifelse(network::is.directed(net), "digraph", "graph")
  if (!is.numeric(arrow.size) || arrow.size < 0) {
    stop("incorrect arrow.size value")
  }
  else if (arrow.size > 0 & is_dir == "graph") {
    warning("network is undirected; arrow.size ignored")
    arrow.size = 0
  }
  if (!is.numeric(arrow.gap) || arrow.gap < 0 || arrow.gap >
      1) {
    stop("incorrect arrow.gap value")
  }
  else if (arrow.gap > 0 & is_dir == "graph") {
    warning("network is undirected; arrow.gap ignored")
    arrow.gap = 0
  }
  if (network::is.hyper(net)) {
    stop("ggnet2 cannot plot hyper graphs")
  }
  if (network::is.multiplex(net)) {
    stop("ggnet2 cannot plot multiplex graphs")
  }
  if (network::has.loops(net)) {
    warning("ggnet2 does not know how to handle self-loops")
  }
  x = max_size
  if (!is.numeric(x) || is.infinite(x) || is.nan(x) || x <
      0) {
    stop("incorrect max_size value")
  }
  data = data.frame(label = get_v(net, "vertex.names"), stringsAsFactors = FALSE)
  data$alpha = set_node(node.alpha, "node.alpha")
  data$color = set_node(node.color, "node.color")
  data$shape = set_node(node.shape, "node.shape")
  data$size = set_node(node.size, "node.size")
  if (length(na.rm) > 1) {
    stop("incorrect na.rm value")
  }
  else if (!is.na(na.rm)) {
    if (!na.rm %in% v_attr) {
      stop(paste("vertex attribute", na.rm, "was not found"))
    }
    x = which(is.na(get_v(net, na.rm)))
    message(paste("na.rm removed", length(x), "nodes out of",
                  nrow(data)))
    if (length(x) > 0) {
      data = data[-x, ]
      network::delete.vertices(net, x)
      if (!nrow(data)) {
        warning("na.rm removed all nodes; nothing left to plot")
        return(invisible(NULL))
      }
    }
  }
  x = size
  if (length(x) == 1 && x %in% c("indegree", "outdegree",
                                 "degree", "freeman")) {

    data$size = sna::degree(net, gmode = is_dir, cmode = ifelse(x == "degree", "freeman", x))

    size.legend = ifelse(is.na(size.legend), x, size.legend)
  }
  x = ifelse(is.na(size.min), 0, size.min)
  if (length(x) > 1 || !is.numeric(x) || is.infinite(x) ||
      is.nan(x) || x < 0) {
    stop("incorrect size.min value")
  }
  else if (x > 0 && !is.numeric(data$size)) {
    warning("node.size is not numeric; size.min ignored")
  }
  else if (x > 0) {
    x = which(data$size < x)
    message(paste("size.min removed", length(x), "nodes out of",
                  nrow(data)))
    if (length(x) > 0) {
      data = data[-x, ]
      network::delete.vertices(net, x)
      if (!nrow(data)) {
        warning("size.min removed all nodes; nothing left to plot")
        return(invisible(NULL))
      }
    }
  }
  x = ifelse(is.na(size.max), 0, size.max)
  if (length(x) > 1 || !is.numeric(x) || is.infinite(x) ||
      is.nan(x) || x < 0) {
    stop("incorrect size.max value")
  }
  else if (x > 0 && !is.numeric(data$size)) {
    warning("node.size is not numeric; size.max ignored")
  }
  else if (x > 0) {
    x = which(data$size > x)
    message(paste("size.max removed", length(x), "nodes out of",
                  nrow(data)))
    if (length(x) > 0) {
      data = data[-x, ]
      network::delete.vertices(net, x)
      if (!nrow(data)) {
        warning("size.max removed all nodes; nothing left to plot")
        return(invisible(NULL))
      }
    }
  }
  x = size.cut
  if (length(x) > 1 || is.null(x) || is.na(x) || is.infinite(x) ||
      is.nan(x)) {
    stop("incorrect size.cut value")
  }
  else if (isTRUE(x)) {
    x = 4
  }
  else if (is.logical(x) && !x) {
    x = 0
  }
  else if (!is.numeric(x)) {
    stop("incorrect size.cut value")
  }
  if (x >= 1 && !is.numeric(data$size)) {
    warning("node.size is not numeric; size.cut ignored")
  }
  else if (x >= 1) {
    x = unique(quantile(data$size, probs = seq(0, 1, by = 1/as.integer(x))))
    if (length(x) > 1) {
      data$size = cut(data$size, unique(x), include.lowest = TRUE)
    }
    else {
      warning("node.size is invariant; size.cut ignored")
    }
  }
  if (!is.null(alpha.palette)) {
    x = alpha.palette
  }
  else if (is.factor(data$alpha)) {
    x = levels(data$alpha)
  }
  else {
    x = unique(data$alpha)
  }
  if (!is.null(names(x))) {
    y = unique(na.omit(data$alpha[!data$alpha %in% names(x)]))
    if (length(y) > 0) {
      stop(paste("no alpha.palette value for", paste0(y,
                                                      collapse = ", ")))
    }
  }
  else if (is.factor(data$alpha) || !is.numeric(x)) {
    data$alpha = factor(data$alpha)
    x = scales::rescale_max(1:length(levels(data$alpha)))
    names(x) = levels(data$alpha)
  }
  alpha.palette = x
  if (!is.null(color.palette)) {
    x = color.palette
  }
  else if (is.factor(data$color)) {
    x = levels(data$color)
  }
  else {
    x = unique(data$color)
  }
  if (length(x) == 1 && "RColorBrewer" %in% rownames(installed.packages()) &&
      x %in% rownames(RColorBrewer::brewer.pal.info)) {
    data$color = factor(data$color)
    n_groups = length(levels(data$color))
    n_colors = RColorBrewer::brewer.pal.info[x, "maxcolors"]
    if (n_groups > n_colors) {
      stop(paste0("too many node groups (", n_groups,
                  ") for ", "ColorBrewer palette ", x, " (max: ",
                  n_colors, ")"))
    }
    else if (n_groups < 3) {
      n_groups = 3
    }
    x = RColorBrewer::brewer.pal(n_groups, x)[1:length(levels(data$color))]
    names(x) = levels(data$color)
  }
  if (!is.null(names(x))) {
    y = unique(na.omit(data$color[!data$color %in% names(x)]))
    if (length(y) > 0) {
      stop(paste("no color.palette value for", paste0(y,
                                                      collapse = ", ")))
    }
  }
  else if (is.factor(data$color) || !is_col(x)) {
    data$color = factor(data$color)
    x = gray.colors(length(x))
    names(x) = levels(data$color)
  }
  color.palette = x
  if (!is.null(shape.palette)) {
    x = shape.palette
  }
  else if (is.factor(data$shape)) {
    x = levels(data$shape)
  }
  else {
    x = unique(data$shape)
  }
  if (!is.null(names(x))) {
    y = unique(na.omit(data$shape[!data$shape %in% names(x)]))
    if (length(y) > 0) {
      stop(paste("no shape.palette value for", paste0(y,
                                                      collapse = ", ")))
    }
  }
  else if (is.factor(data$shape) || !is.numeric(x)) {
    data$shape = factor(data$shape)
    x = (scales::shape_pal())(length(levels(data$shape)))
    names(x) = levels(data$shape)
  }
  shape.palette = x
  if (!is.null(size.palette)) {
    x = size.palette
  }
  else if (is.factor(data$size)) {
    x = levels(data$size)
  }
  else {
    x = unique(data$size)
  }
  if (!is.null(names(x))) {
    y = unique(na.omit(data$size[!data$size %in% names(x)]))
    if (length(y) > 0) {
      stop(paste("no size.palette value for", paste0(y,
                                                     collapse = ", ")))
    }
  }
  else if (is.factor(data$size) || !is.numeric(x)) {
    data$size = factor(data$size)
    x = 1:length(levels(data$size))
    names(x) = levels(data$size)
  }
  size.palette = x
  l = node.label
  if (isTRUE(l)) {
    l = data$label
  }
  else if (length(l) > 1 & length(l) == n_nodes) {
    data$label = l
  }
  else if (length(l) == 1 && l %in% v_attr) {
    l = get_v(net, l)
  }
  else {
    l = ifelse(data$label %in% l, data$label, "")
  }
  if (is.character(mode) && length(mode) == 1) {
    mode = paste0("gplot.layout.", mode)
    if (!exists(mode)) {
      stop(paste("unsupported placement method:", mode))
    }
    xy = network::as.matrix.network.adjacency(net)
    xy = do.call(mode, list(xy, layout.par))
    xy = data.frame(x = xy[, 1], y = xy[, 2])
  }
  else if (is.character(mode) && length(mode) == 2) {
    xy = data.frame(x = set_attr(mode[1]), y = set_attr(mode[2]))
  }
  else if (is.numeric(mode) && is.matrix(mode)) {
    xy = data.frame(x = set_attr(mode[, 1]), y = set_attr(mode[,
                                                               2]))
  }
  else {
    stop("incorrect mode value")
  }
  if (length(mode) == 1) {
    xy$x = scale(xy$x, min(xy$x), diff(range(xy$x)))
    xy$y = scale(xy$y, min(xy$y), diff(range(xy$y)))
  }
  data = cbind(data, xy)
  edges = network::as.matrix.network.edgelist(net)
  if (edge.color[1] == "color" && length(edge.color) == 2) {
    edge.color = ifelse(data$color[edges[, 1]] == data$color[edges[,
                                                                   2]], as.character(data$color[edges[, 1]]), edge.color[2])
    if (!is.null(names(color.palette))) {
      x = which(edge.color %in% names(color.palette))
      edge.color[x] = color.palette[edge.color[x]]
    }
    edge.color[is.na(edge.color)] = edge.color[2]
  }
  edge.color = set_edge(edge.color, "edge.color")
  if (!is_col(edge.color)) {
    stop("incorrect edge.color value")
  }
  edges = data.frame(xy[edges[, 1], ], xy[edges[, 2], ])
  names(edges) = c("X1", "Y1", "X2", "Y2")
  if (!is.null(edge.label)) {
    edges$midX = (edges$X1 + edges$X2)/2
    edges$midY = (edges$Y1 + edges$Y2)/2
    edges$label = set_edge(edge.label, "edge.label")
    edge.label.alpha = set_edge(edge.label.alpha, "edge.label.alpha")
    if (!is.numeric(edge.label.alpha)) {
      stop("incorrect edge.label.alpha value")
    }
    edge.label.color = set_edge(edge.label.color, "edge.label.color")
    if (!is_col(edge.label.color)) {
      stop("incorrect edge.label.color value")
    }
    edge.label.size = set_edge(edge.label.size, "edge.label.size")
    if (!is.numeric(edge.label.size)) {
      stop("incorrect edge.label.size value")
    }
  }
  edge.lty = set_edge(edge.lty, "edge.lty")
  edge.size = set_edge(edge.size, "edge.size")
  if (!is.numeric(edge.size) || any(edge.size <= 0)) {
    stop("incorrect edge.size value")
  }
  p = ggplot(data, aes(x = x, y = y))
  if (nrow(edges) > 0) {
    if (arrow.gap > 0) {
      x.length = with(edges, X2 - X1)
      y.length = with(edges, Y2 - Y1)
      arrow.gap = with(edges, arrow.gap/sqrt(x.length^2 +
                                               y.length^2))
      edges = transform(edges, X1 = X1 + arrow.gap * x.length,
                        Y1 = Y1 + arrow.gap * y.length, X2 = X1 + (1 -
                                                                     arrow.gap) * x.length, Y2 = Y1 + (1 - arrow.gap) *
                          y.length)
    }
    p = p + geom_segment(data = edges, aes(x = X1, y = Y1,
                                           xend = X2, yend = Y2), size = edge.size, color = edge.color,
                         alpha = edge.alpha, lty = edge.lty, arrow = arrow(type = arrow.type,
                                                                           length = unit(arrow.size, "pt")))
  }
  if (nrow(edges) > 0 && !is.null(edge.label)) {
    p = p + geom_point(data = edges, aes(x = midX, y = midY),
                       alpha = edge.alpha, color = edge.label.fill, size = edge.label.size *
                         1.5) + geom_text(data = edges, aes(x = midX,
                                                            y = midY, label = label), alpha = edge.label.alpha,
                                          color = edge.label.color, size = edge.label.size)
  }
  x = list()
  if (is.numeric(data$alpha) && is_one(data$alpha)) {
    x = c(x, alpha = unique(data$alpha))
  }
  if (!is.factor(data$color) && is_one(data$color)) {
    x = c(x, colour = unique(data$color))
  }
  if (is.numeric(data$shape) && is_one(data$shape)) {
    x = c(x, shape = unique(data$shape))
  }
  if (is.numeric(data$size) && is_one(data$size)) {
    x = c(x, size = unique(data$size))
  }
  else {
    x = c(x, size = max_size)
  }
  p = p + geom_point(aes(alpha = factor(alpha), color = factor(color),
                         shape = factor(shape), size = factor(size)))
  if (is.numeric(data$alpha)) {
    v_alpha = unique(data$alpha)
    names(v_alpha) = unique(data$alpha)
    p = p + scale_alpha_manual("", values = v_alpha) + guides(alpha = "none")
  }
  else {
    p = p + scale_alpha_manual(set_name(node.alpha, alpha.legend),
                               values = alpha.palette, breaks = names(alpha.palette),
                               guide = guide_legend(override.aes = x))
  }
  if (!is.null(names(color.palette))) {
    p = p + scale_color_manual(set_name(node.color, color.legend),
                               values = color.palette, breaks = names(color.palette),
                               guide = guide_legend(override.aes = x))
  }
  else {
    v_color = unique(data$color)
    names(v_color) = unique(data$color)
    p = p + scale_color_manual("", values = v_color) + guides(color = "none")
  }
  if (is.numeric(data$shape)) {
    v_shape = unique(data$shape)
    names(v_shape) = unique(data$shape)
    p = p + scale_shape_manual("", values = v_shape) + guides(shape = "none")
  }
  else {
    p = p + scale_shape_manual(set_name(node.shape, shape.legend),
                               values = shape.palette, breaks = names(shape.palette),
                               guide = guide_legend(override.aes = x))
  }
  x = x[names(x) != "size"]
  if (is.numeric(data$size)) {
    v_size = set_size(unique(data$size))
    if (length(v_size) == 1) {
      v_size = as.numeric(names(v_size))
      p = p + scale_size_manual("", values = v_size) +
        guides(size = FALSE)
    }
    else {
      p = p + scale_size_manual(set_name(node.size, size.legend),
                                values = v_size, guide = guide_legend(override.aes = x))
    }
  }
  else {
    p = p + scale_size_manual(set_name(node.size, size.legend),
                              values = set_size(size.palette), guide = guide_legend(override.aes = x))
  }
  if (!is_one(l) || unique(l) != "") {
    label.alpha = set_node(label.alpha, "label.alpha", mode = FALSE)
    if (!is.numeric(label.alpha)) {
      stop("incorrect label.alpha value")
    }
    label.color = set_node(label.color, "label.color", mode = FALSE)
    if (!is_col(label.color)) {
      stop("incorrect label.color value")
    }
    label.size = set_node(label.size, "label.size", mode = FALSE)
    if (!is.numeric(label.size)) {
      stop("incorrect label.size value")
    }
    x = label.trim
    if (length(x) > 1 || (!is.logical(x) & !is.numeric(x) &
                          !is.function(x))) {
      stop("incorrect label.trim value")
    }
    else if (is.numeric(x) && x > 0) {
      l = substr(l, 1, x)
    }
    else if (is.function(x)) {
      l = x(l)
    }
    p = p + geom_text(label = l, alpha = label.alpha, color = label.color,
                      size = label.size, ...)
  }
  x = range(data$x)
  if (!is.numeric(layout.exp) || layout.exp < 0) {
    stop("incorrect layout.exp value")
  }
  else if (layout.exp > 0) {
    x = scales::expand_range(x, layout.exp/2)
  }
  p = p + scale_x_continuous(breaks = NULL, limits = x) +
    scale_y_continuous(breaks = NULL) + theme(panel.background = element_blank(),
                                              panel.grid = element_blank(), axis.title = element_blank(),
                                              legend.key = element_blank(), legend.position = legend.position,
                                              legend.text = element_text(size = legend.size), legend.title = element_text(size = legend.size))
  return(p)
}

