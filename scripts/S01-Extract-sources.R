setwd("~/Shared/Data-Science/Data-Source-Model-Repository/DO/scripts/")

library(git2r)
library(RJSONIO)

desc <- readJSONStream("../DESCRIPTION.json")

sourceFiles <- desc$"source files"
urls <- unlist(lapply(
  sourceFiles,
  function(sf){
    toRet <- sf$"URL template"
    names(toRet) <- sf$"name"
    return(toRet)
  }
))
srcDir <- "../sources/DO"
gitRepo <- urls[1]

## Clone or pull git repository
if(!dir.exists(srcDir)){
  gitRepo <- git2r::clone(url = gitRepo, local_path = srcDir)
}else{
  gitRepo <- git2r::repository(srcDir)
  git2r::pull(gitRepo)
}

###############################################
## Information source files
rcurrent <- git2r::odb_blobs(gitRepo)
rcurrent <- tail(rcurrent[rcurrent$name == "doid.json",], n = 1L)

DO_sourceFiles <- data.frame(url = urls,
                             current = rcurrent$when)

###############################################
## Writing files
toSave <- grep("^DO[_]", ls(), value = T)
ddir <- "../data"

write.table(get(toSave), row.names = FALSE, sep = "\t", quote = FALSE, file=file.path(ddir, paste(toSave, ".txt", sep="")))


