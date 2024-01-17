#' Create the files for labelling the linked clusters from ClusTCR_list_to_matrix function
#' @param my_file Matrix file produce from [ClusTCR]
#' @param max.iter Number of iterations to find the steady state of MCL.
#' @param inflation numeric value
#' @param expansion numeric value
#' @param j part of MCL clustering
#' @import DescTools
#' @export

mcl_cluster <- function(my_file, max.iter=10, inflation = 1, expansion = 1) {
  stat_process <- Sys.time()

  adj.norm <- my_file
  diag(adj.norm) <- 1
  a <- 1
  repeat {
    expans <- adj.norm %^% expansion
    infl <- expans^inflation
    infl.norm <- apply(infl[, ], MARGIN = 2, FUN = function(Var) {
      Var/sum(Var)
    })
    if (identical(infl.norm, adj.norm)) {
      ident <- TRUE
      break
    }
    if (a == max.iter) {
      ident <- FALSE
      a <- a + 1
      break
    }
    adj.norm <- infl.norm
    a <- a + 1
    message(paste("Iteration complete:", a))
  }

  # if (!is.na(infl.norm[1, 1]) & ident) {
  count <- 0
  for (i in 1:ncol(infl.norm)) {
    if (sum(abs(infl.norm[i, ])) != 0) {
      count <- count + 1
    }
    #
  }

  message(paste("Inflation complete"))

  neu <- matrix(nrow = count, ncol = ncol(infl.norm))
  # View(neu)
  zeile <- 1
  for (i in 1:nrow(infl.norm)) {
    if (sum(infl.norm[i, ]) != 0) {
      for (j in 1:ncol(infl.norm)) {
        neu[zeile, j] <- infl.norm[i, j]
      }
      zeile <- zeile + 1
    }
    # message(paste("Processing:", zeile))
  }

  message(paste("MCL complete"))

# Changes neu matrix to 1 and 0
  for (j in 1:ncol(neu)) {
    for (i in 1:nrow(neu)) {
      neu[i,j] <- ifelse(neu[i,j]> 0,1,0)
    }
  }

  message(paste("Finished correcting matrix to binary"))

  for (i in 1:nrow(neu)) {
    num <- ifelse(neu[,i]>1, which(neu[,i] > 1),
                  ifelse((neu[,i] > 1) & neu[i,] == 1 ,which(neu[,i] > 1),
                         ifelse(neu[,i] == 1,i,0)))
    neu[i,] <- num
    neu[,i] <- num
  }

  message(paste("relabelling nodes for MCL start."))

  for (r in 1:6) {
    r = r
    start_time <- Sys.time()
    for (i in 1:nrow(neu)) {
      num_mat <- unique(neu[i,])
      num_mat <- num_mat[order(num_mat)]
      num_mat2 <- num_mat[num_mat>0][1]
      neu[i,] <- ifelse(neu[i,]>0,num_mat2,0)
      neu[,i] <- ifelse(neu[,i]>0,num_mat2,0)
    }
    end_time <- Sys.time()
    diff <- end_time - start_time
    message(paste0("iteration ",r," completed in: ",diff))
  }

  df <- matrix(ncol=3,nrow=nrow(neu))
  for (i in 1:dim(df)[1] ) {
    val <-paste0(which(neu[i,] >0))
    val_name <- paste(val[1])
    val_name
    val_len <- length(paste0(which(neu[i,] >0)))
    for (j in 2:length(val)) {
      val_name <- paste(val_name,val[j])
    }

    num_mat <- unique(neu[i,])
    df[i,1] <- num_mat[num_mat>0][1]
    df[i,2] <- val_name
    df[i,3] <- val_len
  }

  df <- as.data.frame(df)
  df$order <- 1:dim(df)[1]

  # add in
  df2 <- df[c("V1")]
  df2$count <- 1
  df3 <- as.data.frame(ddply(df2,"V1",numcolwise(sum)))
  df3 <- df3[order(as.numeric(df3$count),decreasing = T), ]
  df3$Clust_size_order <- 1:dim(df3)[1]

  df_num_rep <- as.data.frame(unique(df[,1]))
  names(df_num_rep) <- "V1"
  df_num_rep$cluster <- 1:dim(df_num_rep)[1]

  df <- merge(df,df_num_rep,by="V1",sort=F)
  df <- merge(df,df3, by=c("V1"))
  # names(df)[1:3] <- c("original_cluster","Cluster size","Nodes")

  df <- df[order(df$order),]
  names(df)[1:3] <- c("Original_cluster","nodes","#_of_connections")
  df$CDR3_Vgene <- colnames(my_file)

  stat_end_process <- Sys.time()
  diff <- stat_end_process - stat_process
  message(paste("this process took", round(diff,2),"seconds"))
  message(paste("Completed process, can now display data as either motif or netplot."))
  mylist<-list(Cluster_lab = df,
               Normalised_tabel = infl.norm)
}


#' Create the files for labelling the linked clusters from ClusTCR_list_to_matrix function
#' @param my_file Matrix file produce from [ClusTCR]
#' @param max.iter Number of iterations to find the steady state of MCL.
#' @param inflation numeric value
#' @param expansion numeric value
#' @param j part of MCL clustering
#' @import DescTools
#' @export

mcl_cluster_large <- function(my_file, max.iter=10, inflation = 1, expansion = 1) {
  df_mat2 <- my_file

  num_large_issue <- length(df_mat2)
  Cluster_df <- as.list(NULL)
  infl.norm_mat <- as.list(NULL)
  for (k in 1:num_large_issue) {
    message(paste("Starting cluster",k,"of",num_large_issue))
    # message(paste("Starting",k,"of",num))
    df_mat3 <- df_mat2[[k]]
    adj.norm <- df_mat3
    diag(adj.norm) <- 1
    a <- 1
    repeat {
      expans <- adj.norm %^% expansion
      infl <- expans^inflation
      infl.norm <- apply(infl[, ], MARGIN = 2, FUN = function(Var) {
        Var/sum(Var)
      })
      if (identical(infl.norm, adj.norm)) {
        ident <- TRUE
        break
      }
      if (a == max.iter) {
        ident <- FALSE
        a <- a + 1
        break
      }
      adj.norm <- infl.norm
      a <- a + 1
    }

    count <- 0
    for (i in 1:ncol(infl.norm)) {
      if (sum(abs(infl.norm[i, ])) != 0) {
        count <- count + 1
      }
      #
    }

    message(paste("Inflation complete for",k))

    neu <- matrix(nrow = count, ncol = ncol(infl.norm))
    # View(neu)
    zeile <- 1
    for (i in 1:nrow(infl.norm)) {
      if (sum(infl.norm[i, ]) != 0) {
        for (j in 1:ncol(infl.norm)) {
          neu[zeile, j] <- infl.norm[i, j]
        }
        zeile <- zeile + 1
      }
      # message(paste("Processing:", zeile))
    }

    message(paste("MCL complete for",k))

    # Changes neu matrix to 1 and 0
    for (j in 1:ncol(neu)) {
      for (i in 1:nrow(neu)) {
        neu[i,j] <- ifelse(neu[i,j]> 0,1,0)
      }
    }

    message(paste("Finished correcting matrix to binary"))

    for (i in 1:nrow(neu)) {
      num <- ifelse(neu[,i]>1, which(neu[,i] > 1),
                    ifelse((neu[,i] > 1) & neu[i,] == 1 ,which(neu[,i] > 1),
                           ifelse(neu[,i] == 1,i,0)))
      neu[i,] <- num
      neu[,i] <- num
    }

    message(paste("relabelling nodes for MCL start.",k))

    for (r in 1:6) {
      r = r
      start_time <- Sys.time()
      for (i in 1:nrow(neu)) {
        num_mat <- unique(neu[i,])
        num_mat <- num_mat[order(num_mat)]
        num_mat2 <- num_mat[num_mat>0][1]
        neu[i,] <- ifelse(neu[i,]>0,num_mat2,0)
        neu[,i] <- ifelse(neu[,i]>0,num_mat2,0)
      }
      end_time <- Sys.time()
      diff <- end_time - start_time
      # message(paste0("iteration ",r," completed in: ",diff))
    }

    df <- matrix(ncol=3,nrow=nrow(neu))

    for (i in 1:dim(df)[1] ) {
      val <-paste0(which(neu[i,] >0))
      val_name <- paste(val[1])
      val_name
      val_len <- length(paste0(which(neu[i,] >0)))
      for (j in 2:length(val)) {
        val_name <- paste(val_name,val[j])
      }

      num_mat <- unique(neu[i,])
      df[i,1] <- num_mat[num_mat>0][1]
      df[i,2] <- val_name
      df[i,3] <- val_len
    }

    df <- as.data.frame(df)
    df$order <- 1:dim(df)[1]

    # add in
    df2 <- df[c("V1")]
    df2$count <- 1
    df3 <- as.data.frame(ddply(df2,"V1",numcolwise(sum)))
    df3 <- df3[order(as.numeric(df3$count),decreasing = T), ]
    df3$Clust_size_order <- 1:dim(df3)[1]
    df3$list_pos <- k
    df3$Clust_size_order_k <- paste(df3$Clust_size_order,k,sep="_")
    df_num_rep <- as.data.frame(unique(df[,1]))
    names(df_num_rep) <- "V1"
    df_num_rep$cluster <- 1:dim(df_num_rep)[1]

    df <- merge(df,df_num_rep,by="V1",sort=F)
    df <- merge(df,df3, by=c("V1"))

    df <- df[order(df$order),]
    names(df)[1:3] <- c("Original_cluster","nodes","#_of_connections")
    df$Original_cluster <- paste(df$Original_cluster,k,sep="_")
    df$CDR3_Vgene <- colnames(df_mat3)
    Cluster_df[[k]] <- df
    infl.norm_mat[[k]] <- infl.norm
  }

  Cluster_df2 <- Cluster_df[[1]]

  for (i in 2:length(Cluster_df)) {
    Cluster_df2 <- rbind(Cluster_df2,Cluster_df[[i]])
    message("combinded",i,"of",length(Cluster_df))
  }
  dim(Cluster_df2)
  Cluster_df2 <- Cluster_df2[order(Cluster_df2$count,decreasing = T),]

  mylist<-list(Cluster_lab = Cluster_df2,
               Normalised_tabel = infl.norm_mat)
  mylist
}
