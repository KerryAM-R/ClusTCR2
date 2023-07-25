#' Creates ClusTCR matrix
#' This function identifies similar CDR3 amino acid sequences based on the same length and V_gene
#' @param my_file uploaded file with junction_aa (CD3 sequences), variable gene.
#' @param v_gene Variable gene column name
#' @param allele The allele, if present as *00 will be removed if the user requires it.
#' @param cores_selected Chose # of cores for parallel functions
#' @return X by Y matrix of strucutrally related CDR3 sequences.
#' @importFrom stringr str_split_fixed
#' @importFrom doParallel registerDoParallel
#' @import plyr
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @import DescTools
#' @importFrom grDevices rainbow
#' @importFrom stats setNames
#' @importFrom utils head
#' @export

ClusTCR <- function(my_file, allele=NULL, v_gene = "v_call", cores_selected = 4) {

  registerDoParallel(cores = cores_selected)
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
    unique(amino_acid_test_top2$V_call_len)
    df_len <- as.data.frame(unique(amino_acid_test_top2$V_call_len))
    names(df_len) <- "Len"
    df_len <- as.data.frame(df_len[order(df_len$Len),])
    names(df_len) <- "Len"
    df_len
    edge <- as.data.frame(matrix(nrow = 1, ncol = 3))
    names(edge) <- c("source", "target", "Val")
    head(edge)

    message("creating empty matrixes")
    res.all <- foreach(j=1:dim(df_len)[1]) %dopar% {
      df.clust_1 <- subset(amino_acid_test_top2,amino_acid_test_top2$V_call_len==df_len[j,1])
      if (length(df.clust_1)>1) {
        res.all <- as.data.frame(matrix(nrow = dim(df.clust_1)[1], ncol =  dim(df.clust_1)[1]))
        rownames(res.all) <- df.clust_1$Vgene_cdr3
        names(res.all) <- df.clust_1$Vgene_cdr3
        res.all
      }
    }
    message("Performing edit distance")
    sim2 <- foreach(j=1:dim(df_len)[1]) %dopar% {
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
      sim2 <- res.all[[j]]
      message(paste("Completed matrix",j))
    }

    f <- function(m) {
      m[lower.tri(m)] <- t(m)[lower.tri(m)]
      m
    }
      fsim2 <- foreach(j=1:dim(df_len)[1]) %dopar% {
      sim2 <- f(sim2[[j]])
    }

      message(paste("keeping edit distance of 1"))
    ham.vals <- foreach(j=1:dim(df_len)[1]) %dopar% {
      ham.vals <- setNames(
        cbind(
          rev(expand.grid(rownames(fsim2[[j]]), names(fsim2[[j]]))),
          c(t(fsim2[[j]]))
        ), c("source", "target", "Val")
      )
      ham.vals_2 <- subset(ham.vals,ham.vals$Val==1)
    }

    message(paste("Creating target and source object"))

    df_net3 <- as.data.frame((ham.vals[[1]][1:2]))
    head(df_net3)
    for (i in 2:dim(df_len)[1]) {
      df_net2 <- as.data.frame((ham.vals[[i]][1:2]))
      df_net3 <- rbind(df_net3,df_net2)
    }

    df_net3$count <- 1
    df_net3 <- df_net3[order(df_net3$source),]
    message(paste("Creating matrix for MCL"))
    df_mat <- (table(as.character(df_net3$source), as.character(df_net3$target)))
    message(paste("Matrix complete"))
    df_mat
  }
  else {
    message("Incorrect V gene column")
  }
}
