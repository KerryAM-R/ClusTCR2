#' Creates ClusTCR matrix
#' This function identifies similar CDR3 amino acid sequences based on the same length and V_gene
#' @param my_file uploaded file with junction_aa (CD3 sequences), variable gene.
#' @param v_gene Variable gene column name
#' @param allele The allele, if present as *00 will be removed if the user requires it.
#' @return X by Y matrix of strucutrally related CDR3 sequences.
#' @importFrom stringr str_split_fixed
#' @importFrom plyr ddply
#' @importFrom DescTools StrDist
#' @importFrom stats setNames
#' @importFrom utils head
#' @export

ClusTCR <- function(my_file, allele=NULL, v_gene = "v_call") {
  amino_acid_test_top <- my_file
  amino_acid_test_top2 <- amino_acid_test_top[!duplicated(amino_acid_test_top$junction_aa), ]
  amino_acid_test_top2$len <- nchar(amino_acid_test_top2$junction_aa)

  if (is.null(allele)) {
    stop("allele has to be TRUE or FALSE")
  }

  if (length(amino_acid_test_top2[,names(amino_acid_test_top2) %in% "v_call"]) > 0) {
    amino_acid_test_top2$v_call <- gsub("[*]..","",amino_acid_test_top2[,names(amino_acid_test_top2) %in% v_gene])
    # v_call <- as.data.frame(str_split_fixed(amino_acid_test_top2$v_call, '-', 2))
    v_call <- as.data.frame(amino_acid_test_top2$v_call)
    names(v_call) <- "V1"
    amino_acid_test_top2$v_call <- v_call$V1
    amino_acid_test_top2$count <- 1
    amino_acid_test_top2$Vgene_cdr3 <- paste(amino_acid_test_top2$junction_aa,amino_acid_test_top2$v_call,sep = "_")
    amino_acid_test_top2$V_call_len <- paste(amino_acid_test_top2$len,amino_acid_test_top2$v_call,sep = "_")

    df_len <- as.data.frame(unique(amino_acid_test_top2$V_call_len))
    names(df_len) <- "Len"
    df_len <- as.data.frame(df_len[order(df_len$Len),])
    names(df_len) <- "Len"

    message("creating empty matrixes")
    res.all <- as.list(NULL)
    for(j in 1:dim(df_len)[1]) {
      df.clust_1 <- subset(amino_acid_test_top2,amino_acid_test_top2$V_call_len==df_len[j,1])
      # print(dim(df.clust_1)[1])

      if (dim(df.clust_1)[1]>1) {
        res.all[[j]] <- as.data.frame(matrix(nrow = dim(df.clust_1)[1], ncol =  dim(df.clust_1)[1]))
        rownames(res.all[[j]]) <- df.clust_1$Vgene_cdr3
        names(res.all[[j]]) <- df.clust_1$Vgene_cdr3
        res.all
      }
    }

    res.all <- res.all[!sapply(res.all,is.null)]
    res.all
    sim2 <- as.list(NULL)
    message("Performing edit distance")
    length(res.all)
    for (j in 1:length(res.all))  {
      df <- as.data.frame(res.all[[j]])
      for(r in 1:dim(df)[1]) {
        for (i in 1:dim(df)[1]) {
          if (r == i) {
          }
          else if (i>r) {
          }
          else {
            res <- StrDist(rownames(df)[i], names(df)[r], method = "hamming", mismatch = 1, gap = 1, ignore.case = FALSE)
            res.all[[j]][i,r] <- as.numeric(res)
          }
        }
      }
      sim2[[j]] <- res.all[[j]]
    }

    f <- function(m) {
      m[lower.tri(m)] <- t(m)[lower.tri(m)]
      m
    }

    length(res.all)
    fsim2 <- as.list(NULL)

    for (j in 1:length(res.all)) {
      fsim2[[j]] <- f(sim2[[j]])
    }

    fsim2

    message(paste("keeping edit distance of 1"))


    ham.vals <- as.list(NULL)

    for (j in 1:length(fsim2)) {
      ham.vals2 <- setNames(
        cbind(
          rev(expand.grid(rownames(fsim2[[j]]), names(fsim2[[j]]))),
          c(t(fsim2[[j]]))
        ), c("source", "target", "Val")
      )

      setNames(
        cbind(
          rev(expand.grid(rownames(fsim2[[j]]), names(fsim2[[j]]))),
          c(t(fsim2[[j]]))
        ), c("source", "target", "Val")
      )
      ham.vals[[j]] <-  ham.vals2
    }

    ham.vals2 <- as.list(NULL)

    for (i in 1:length(ham.vals)) {

      if (dim(subset(ham.vals[[i]],ham.vals[[i]]$Val==1))[1]>0) {
        ham.vals2[[i]] <- subset(ham.vals[[i]],ham.vals[[i]]$Val==1)
      }
      else {
        ham.vals2[[i]] <- NULL

      }
    }

    ham.vals2 <- ham.vals2[!sapply(ham.vals2,is.null)]

    if ( length(ham.vals2) == 0 ) {
      message("No clusters == 1 edit distance and therefore MCL not performed")
      as.data.frame("No clusters found")
    }

    else {
      message(paste("Creating target and source object"))

      df_net3 <- as.data.frame((ham.vals2[[1]][1:2]))

      for (i in 2:length(ham.vals2)[1]) {
        df_net2 <- as.data.frame((ham.vals2[[i]][1:2]))
        df_net3 <- rbind(df_net3,df_net2)
      }

      df_net3$count <- 1
      df_net3 <- df_net3[order(df_net3$source),]
      message(paste("Creating matrix for MCL"))
      df_mat <- (table(as.character(df_net3$source), as.character(df_net3$target)))
      message(paste("Matrix complete"))
      df_mat
    }
  }
  else {
    message("Incorrect V gene column")
  }
}

#' Creates ClusTCR matrix
#' This function identifies similar CDR3 amino acid sequences based on the same length and V_gene
#' @param my_file uploaded file with junction_aa (CD3 sequences), variable gene.
#' @param v_gene Variable gene column name
#' @param allele The allele, if present as *00 will be removed if the user requires it.
#' @return X by Y matrix of strucutrally related CDR3 sequences.
#' @importFrom stringr str_split_fixed
#' @importFrom DescTools StrDist
#' @importFrom plyr ddply
#' @importFrom stats setNames
#' @importFrom utils head
#' @export

ClusTCR_Large <- function(my_file, allele=NULL, v_gene = "v_call") {
  amino_acid_test_top <- my_file
  amino_acid_test_top2 <- amino_acid_test_top[!duplicated(amino_acid_test_top$junction_aa), ]
  amino_acid_test_top2$len <- nchar(amino_acid_test_top2$junction_aa)

  if (is.null(allele)) {
    stop("allele has to be TRUE or FALSE")
  }

  if (length(amino_acid_test_top2[,names(amino_acid_test_top2) %in% "v_call"]) > 0) {
    amino_acid_test_top2$v_call <- gsub("[*]..","",amino_acid_test_top2[,names(amino_acid_test_top2) %in% v_gene])
    # v_call <- as.data.frame(str_split_fixed(amino_acid_test_top2$v_call, '-', 2))
    v_call <- as.data.frame(amino_acid_test_top2$v_call)
    names(v_call) <- "V1"
    amino_acid_test_top2$v_call <- v_call$V1
    amino_acid_test_top2$count <- 1
    amino_acid_test_top2$Vgene_cdr3 <- paste(amino_acid_test_top2$junction_aa,amino_acid_test_top2$v_call,sep = "_")
    amino_acid_test_top2$V_call_len <- paste(amino_acid_test_top2$len,amino_acid_test_top2$v_call,sep = "_")

    df_len <- as.data.frame(unique(amino_acid_test_top2$V_call_len))
    names(df_len) <- "Len"
    df_len <- as.data.frame(df_len[order(df_len$Len),])
    names(df_len) <- "Len"

    message("creating empty matrixes")
    res.all <- as.list(NULL)
    for(j in 1:dim(df_len)[1]) {
      df.clust_1 <- subset(amino_acid_test_top2,amino_acid_test_top2$V_call_len==df_len[j,1])
      if (dim(df.clust_1)[1]>1) {
        res.all[[j]] <- as.data.frame(matrix(nrow = dim(df.clust_1)[1], ncol =  dim(df.clust_1)[1]))
        rownames(res.all[[j]]) <- df.clust_1$Vgene_cdr3
        names(res.all[[j]]) <- df.clust_1$Vgene_cdr3
        res.all
      }
    }

    res.all <- res.all[!sapply(res.all,is.null)]
    sim2 <- as.list(NULL)

    length(res.all)
    for (j in 1:length(res.all))  {
      message(paste("starting edit distance on",j,"that has",dim(res.all[[j]])[2],"unique sequences"))
      df <- as.data.frame(res.all[[j]])
      for(r in 1:dim(df)[1]) {
        for (i in 1:dim(df)[1]) {
          if (r == i) {
          }
          else if (i>r) {
          }
          else {
            res <- StrDist(rownames(df)[i], names(df)[r], method = "hamming", mismatch = 1, gap = 1, ignore.case = FALSE)
            res.all[[j]][i,r] <- as.numeric(res)
          }
        }
      }
      sim2[[j]] <- res.all[[j]]
    }

    f <- function(m) {
      m[lower.tri(m)] <- t(m)[lower.tri(m)]
      m
    }

    length(res.all)
    fsim2 <- as.list(NULL)

    for (j in 1:length(res.all)) {
      fsim2[[j]] <- f(sim2[[j]])
    }
    message(paste("keeping edit distance of 1"))
    ham.vals <- as.list(NULL)

    for (j in 1:length(fsim2)) {
      message(paste("Creating source-target files",j,"of",length(fsim2)))
      ham.vals2 <- setNames(
        cbind(
          rev(expand.grid(rownames(fsim2[[j]]), names(fsim2[[j]]))),
          c(t(fsim2[[j]]))
        ), c("source", "target", "Val")
      )

      setNames(
        cbind(
          rev(expand.grid(rownames(fsim2[[j]]), names(fsim2[[j]]))),
          c(t(fsim2[[j]]))
        ), c("source", "target", "Val")
      )
      ham.vals[[j]] <-  ham.vals2
    }

    ham.vals2 <- as.list(NULL)

    for (i in 1:length(ham.vals)) {

      if (dim(subset(ham.vals[[i]],ham.vals[[i]]$Val==1))[1]>0) {
        ham.vals2[[i]] <- subset(ham.vals[[i]],ham.vals[[i]]$Val==1)
      }
      else {
        ham.vals2[[i]] <- NULL

      }
    }

    ham.vals2 <- ham.vals2[!sapply(ham.vals2,is.null)]

    if ( length(ham.vals2) == 0 ) {
      message("No clusters == 1 edit distance and therefore MCL not performed")
      as.data.frame("No clusters found")
    } else {
      df_mat2 <- as.list(NULL)
      num.ham.vals2 <- length(ham.vals2)

      message(paste("Creating target and source object"))
      for (i in 1:num.ham.vals2) {
        df_net3 <- as.data.frame((ham.vals2[[i]][1:2]))
        df_net3$count <- 1
        df_net3 <- df_net3[order(df_net3$source),]
        message(paste("Creating matrix for MCL",i))
        df_mat <- (table(as.character(df_net3$source), as.character(df_net3$target)))
        df_mat2[[i]] <- df_mat
      }

      df_mat2

    }
  } else {
    message("Incorrect V gene column")
  }
}


#' @example example.R
