
SNPdistance <-function(input, missing.threshold = 80, includeHet=FALSE){
  if (!require(seqinr)){
    install.packages("seqinr",repos = "http://cran.us.r-project.org")
    library(seqinr)
  }
  options(stringsAsFactors = F)
  
  ## matrix from FASTA
  dna <- read.fasta(file = input,forceDNAtolower = F)
  dna_matrix<-t(matrix(unlist(dna),ncol = length(dna)))
  rownames(dna_matrix)<-names(dna)
  if (!includeHet){
    dna_matrix[which(dna_matrix!="A" & dna_matrix!="C" & dna_matrix!="G" & dna_matrix!="T")]<-NA
  } else {
    dna_matrix[which(dna_matrix=="N")]<-NA
  }
  
  ## distance
  SNPdel <- apply(dna_matrix, 2, function(x) ((length(which(is.na(x))==T)/nrow(dna_matrix))*100)>missing.threshold)
  dna_matrix <- dna_matrix[, !SNPdel]
  res <- matrix(0,ncol=3,nrow = (nrow(dna_matrix) * (nrow(dna_matrix) - 1)/2))
  num <- 1L
  for (i in 1:(nrow(dna_matrix) - 1)) {
    for (j in (i+1):nrow(dna_matrix)) {
      com<- dna_matrix[i, ] != dna_matrix[j, ]
      res[num,1]<-row.names(dna_matrix)[i]
      res[num,2]<-row.names(dna_matrix)[j]
      res[num,3] <-  sum(com, na.rm = TRUE)
      num <- num + 1L
    }
  }
  res
  write.csv(res,paste0(unlist(strsplit(input, split='.fasta', fixed=TRUE))[1],"_SNPdistance.csv"),row.names = F)
}

