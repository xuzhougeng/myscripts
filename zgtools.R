# fix the dependent problem
fixit <- function(){
  
  # set the language to english
  env <- Sys.getenv("LANGUAGE")
  on.exit(Sys.setenv("LANGUAGE" = env))
  Sys.setenv("LANGUAGE" = "en")
  
  # get the history
  file1 <- tempfile("Rrawhist")
  savehistory(file1)
  rawhist <- readLines(file1)
  unlink(file1)

  # get the package name
  rawhist <- rawhist[grepl("^library|^require", rawhist)]
  rawhist <- rawhist[length(rawhist)]
  package.name <- strsplit(rawhist, "\\(|\\)")[[1]][2]
  package.name <- gsub("'|\"","", package.name)
  
  # get the dependent package
  text <- try(library(package.name, character.only = TRUE))
  
  # extract the package name  from error message
  package.name <- gsub(".*called ‘(.*?)’.*", "\\1", text[1])
  
  # download package
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(package.name)
  
}

# read the STAR count result and build a matrix
read.STAR <- function(dir){
	  files <- list.files(dir, pattern = "*_ReadsPerGene.out.tab", full.names = TRUE)
  
  file <- names <- gsub("_ReadsPerGene.out.tab","",basename(files))
    file <- names <- gsub("-"," <- ", file <- names)
    
    df.list <- lapply(files, function(x) {
						      data.table::fread(x, skip=4)[, c("V1", "V2")]
							    }
	  )
	  df <- do.call(cbind,lapply(df.list, function(x){
									     x[,2]
										   }
	  ))
	  mt <- as.matrix(df)
	    row.names(mt) <- df.list[[1]]$V1
	    colnames(mt) <- file <- names
		  mt
		  
}

# Plot PCA without building a DESeq2 object
plotPCA <- function(mt,  ntop = 500, group = NULL, blind = FALSE) {
  
  vsd <- DESeq2::vst(mt, blind = blind)
  
  rv <- matrixStats::rowVars( vsd )
  
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(vsd[select, ]))
  
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  if (is.null(group)){
    
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
    p <- ggplot(data = d, aes_string(x = "PC1", y = "PC2")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] *  100), "% variance")) + 
      coord_fixed()
    return(p)
  } else{
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group)
    p <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color= "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] *  100), "% variance")) + 
      coord_fixed()
    return(p)
  }

}
