#' Create the files for labelling the linked clusters from ClusTCR_list_to_matrix function
#' @param my_file Matrix file produce from [ClusTCR_list_to_matrix]
#' @param max.iter Number of iterations to find the steady state of MCL.
#' @param inflation numeric value
#' @param expansion numeric value
#' @param j part of MCL clustering
#' @import DescTools
#' @export


mcl_cluster <- function(my_file, max.iter=100, inflation = 1, expansion = 1) {
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
  }



  # if (!is.na(infl.norm[1, 1]) & ident) {
  count <- 0
  for (i in 1:ncol(infl.norm)) {
    if (sum(abs(infl.norm[i, ])) != 0) {
      count <- count + 1
    }
  }

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
  }

# Changes neu matrix to 1 and 0
  for (j in 1:ncol(neu)) {
    for (i in 1:nrow(neu)) {
      neu[i,j] <- ifelse(neu[i,j]> 0,1,0)
    }
  }
  nrow(neu)
  ncol(neu)

  for (i in 1:nrow(neu)) {
    num <- ifelse(neu[,i]>1, which(neu[,i] > 1),
                  ifelse((neu[,i] > 1) & neu[i,] == 1 ,which(neu[,i] > 1),
                         ifelse(neu[,i] == 1,i,0)))
    neu[i,] <- num
    neu[,i] <- num
  }

  for (i in 1:nrow(neu)) {
    num_mat <- unique(neu[i,])
    num_mat <- num_mat[order(num_mat)]
    num_mat2 <- num_mat[num_mat>0][1]
    neu[i,] <- ifelse(neu[i,]>0,num_mat2,0)
    neu[,i] <- ifelse(neu[,i]>0,num_mat2,0)
  }

  for (i in 1:nrow(neu)) {
    num_mat <- unique(neu[i,])
    num_mat <- num_mat[order(num_mat)]
    num_mat2 <- num_mat[num_mat>0][1]
    neu[i,] <- ifelse(neu[i,]>0,num_mat2,0)
    neu[,i] <- ifelse(neu[,i]>0,num_mat2,0)
  }

  for (i in 1:nrow(neu)) {
    num_mat <- unique(neu[i,])
    num_mat <- num_mat[order(num_mat)]
    num_mat2 <- num_mat[num_mat>0][1]
    neu[i,] <- ifelse(neu[i,]>0,num_mat2,0)
    neu[,i] <- ifelse(neu[,i]>0,num_mat2,0)
  }

  for (i in 1:nrow(neu)) {
    num_mat <- unique(neu[i,])
    num_mat <- num_mat[order(num_mat)]
    num_mat2 <- num_mat[num_mat>0][1]
    neu[i,] <- ifelse(neu[i,]>0,num_mat2,0)
    neu[,i] <- ifelse(neu[,i]>0,num_mat2,0)
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
  mylist<-list(Cluster_lab = df,
               Normalised_tabel = infl.norm)


}
