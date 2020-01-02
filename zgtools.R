# fix the dependent problem
fix <- function(){
  
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
