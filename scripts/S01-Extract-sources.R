library(git2r)
library(RJSONIO)
library(here)
source(here("../00-Utils/downloadSourceFiles.R"))

# reads json file "Description"
desc <- readJSONStream(here("DESCRIPTION.json"))

# extracts source file name (.owl) and url from Json Description file
sourceFiles <- desc$"source files" 
urls <- unlist(lapply(
  sourceFiles,
  function(sf){
    toRet <- sf$"URL template"
    names(toRet) <- sf$"name"
    return(toRet)
  }
))
srcDir <- here("sources")

downloadSourceFiles(urls, srcDir, httpForce = TRUE)

# ## Clone or pull git repository 
# srcDir <- here("sources/HumanDiseaseOntology") # source directory
# gitRepo <- urls[1] # not useful for updates, only for first pull
# if(!dir.exists(srcDir)){ # first pull
#   gitRepo <- clone(url = gitRepo, local_path = srcDir)
# }else{ # updates
#   gitRepo <- repository(srcDir) # opens source directory
#   fetch(repo = gitRepo, name = "origin") # pulls new data --> new data are in source directory as well
#   checkout(gitRepo, "master")
# }
# 
# # Check if some files in sources>DO>src>ontology have been updated
# 
# ###############################################
# ## Extract from new pull the time and date --> saves on a new D0_sourceFiles, with url from json Description file
# rcurrent <- git2r::odb_blobs(gitRepo)
# rcurrent <- tail(rcurrent[rcurrent$name == "doid.json",], n = 1L)
# DO_sourceFiles <- data.frame(url = urls,
#                              current = rcurrent$when)
# 
# # lf <- grep("doid.json", list.files(file.path(srcDir, "src", "ontology"), full.names = T), value = T)
# # DO_sourceFiles <- data.frame(url = urls,
# #                              current = file.info(lf)$mtime)
# 
# ###############################################
# ## Saves new D0_sourceFiles.txt (only contains url and time/date of pull) in Data
# toSave <- grep("^DO[_]", ls(), value = T) # selects D0_sourceFiles
# ddir <- "../data" # defines where to save new data
# write.table(get(toSave), row.names = FALSE, sep = "\t", quote = FALSE, file=file.path(ddir, paste(toSave, ".txt", sep="")))


