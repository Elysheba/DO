rm(list = ls())
gc()


library(XML)
library(parallel)
library(git2r)
library(RJSONIO)
# source(here("../00-Utils/writeLastUpdate.R"))
library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(ReDaMoR)
##
mc.cores <- 20
# sdir <- here("sources")
ddir <- here("data")

###############################################################################@
## Source information - extract source file name (.owl) from description JSON file
###############################################################################@

desc <- RJSONIO::readJSONStream(here("DESCRIPTION.json"))
sourceFiles <- desc$"source files"
sfi_name <- unlist(lapply(
  sourceFiles,
  function(sf){
    toRet <- sf$"name"
    return(toRet)
  }
))
sdir <- here("sources", gsub(".zip", "", sfi_name), "src/ontology")
# sdir <- here("sources")

###############################################################################@
## Data from do ----
###############################################################################@
## decompress gz
# for(f in sfi_name$file){
  gzf <- file.path(sdir,sfi_name)
  print(gzf)
  dest <- gsub("-20.*", "", gzf)
  system(paste0("unzip ", gzf, " -d ", sdir))
  # system(paste0("rm ",file.path(sdir,gsub(".gz","",f))))
# }


###############################################################################@
## Data model
###############################################################################@
# load(here("model", "DO.rda"))
# dm <- model_relational_data()
# save(dm, file = here("model", "DO.rda"))

###############################################################################@
## Data from doid.json ----
###############################################################################@
## Convert OWL to JSON
# uses a program that is not in RStudio but should be run on the shell --> adds in the "shell" path where the program is
# Sys.setenv(PATH = paste(Sys.getenv("PATH"),"~/Shared/Data-Science/Data-Source-Model-Repository/00-Utils/bin/",sep = ":"))
# # Then runs the program to convert the owl source file into a json
# system(paste("robot convert --input ",file.path(sdir,"doid.owl"),
#              " --output ",file.path(sdir,"doid.json"), sep = ""))
# then reads the json filr
readJson <- jsonlite::fromJSON(txt = file.path(sdir,"doid.json"))


###########################################
## from json file, extract nodes (= diseases), with the following: id, def (longer definition), label (name), xref (cross-references in other DBs), syn (synonyms)
nodesJson <- lapply(1:nrow(readJson$graphs$nodes[[1]]),
                    function(i){
                      # print(i)
                      ## id
                      id <- gsub("_",":",gsub(".*/","",readJson$graphs$nodes[[1]]$id[[i]]))
                      ## def
                      def <- readJson$graphs$nodes[[1]]$meta$definition[i, "val"]
                      ## Xref
                      Xref <- readJson$graphs$nodes[[1]]$meta$xrefs[[i]]$val
                      ## name
                      name <- readJson$graphs$nodes[[1]]$meta$synonyms[[i]]$val
                      ## lbl
                      lbl <- readJson$graphs$nodes[[1]]$lbl[[i]]
                      ## alt id
                      if(!is.null(readJson$graphs$nodes[[1]]$meta$basicPropertyValues[[i]])){
                        alt <- readJson$graphs$nodes[[1]]$meta$basicPropertyValues[[i]] %>%
                          filter(grepl("hasAlternativeId", pred)) %>%
                          pull(val)
                      }else{
                        alt <- NULL
                      }
                      ## df
                      df1 <- tibble(
                        id = id,
                        def = def,
                        label = lbl)
                      if(length(Xref) == 0){
                        df2 <- NULL
                      }else{
                        df2 <- tibble(
                          id = id,
                          xref = Xref)
                      }
                      if(length(name) == 0){
                        df3 <- NULL
                      }else{
                        df3 <- tibble(
                          id = id,
                          syn = name)
                      }
                      if(length(alt) == 0){
                        df4 <- NULL
                      }else{
                        df4 <- tibble(
                          id = id, 
                          alt = alt
                        )
                      }
                      return(list(id = df1, xref = df2, syn = df3, alt = df4))
                    }
)

# Generate tables with id, cross-references and synonyms only
id <- do.call(rbind, lapply(nodesJson, function(x) x$id))
xref <- do.call(rbind, lapply(nodesJson, function(x) x$xref))
syn <- do.call(rbind, lapply(nodesJson, function(x) x$syn))
alt <- do.call(rbind, lapply(nodesJson, function(x) x$alt))

###########################################
## extract edges (parents: ... "is a " ...) 
edgesJson <- readJson$graphs$edges[[1]]
edgesJson <- edgesJson[which(edgesJson$pred %in% c("is_a")),]
# simplify names of subjects and objects
edgesJson$sub <- gsub("_",":",gsub(".*/","",edgesJson$sub))
edgesJson$obj <- gsub("_",":",gsub(".*/","",edgesJson$obj))

# function to get all descendants from a single entity in the tree
getDescendants <- function(sp){
  direct <- edgesJson[which(edgesJson$obj == sp),"sub"]
  descendants <- direct
  level <- 0
  dLev <- c()
  for(d in direct){
    dDesc <- getDescendants(d)
    dLev <- c(dLev, dDesc$level)
    descendants <- c(descendants, dDesc$descendants)
  }
  if(length(dLev)>0){
    level <- max(dLev)+1
  }
  return(list(descendants=unique(c(descendants,sp)), level=level))
}
# get all descendants from DOID:4, which is "disease" (general concept), because DO also contains ontologies for other things than diseases (eg phenotype)
disease <- getDescendants("DOID:4")
#lapply(disease,length)
#dim(edgesJson)
# check
"DOID:4" %in% disease$descendants
# length of disease ids
length(disease$descendants) 

######################################
## crossId 
#dim(xref)
# Only select crossIds from diseases
crossId <- xref[xref$id %in% disease$descendants,]
#dim(crossId)
names(crossId) <- c("dbid1","dbid2")

## Remove spaces
head(grep(": ",crossId$dbid2,value = T))
head(grep(": ",crossId$dbid1,value = T))
crossId$dbid1 <- gsub(" ","", crossId$dbid1)
crossId$dbid2 <- gsub(" ","", crossId$dbid2)
dim(crossId)

# check no issue in names
grep("#",crossId$dbid1,value = T)
grep("#",crossId$dbid2,value = T)
# table with number of different databases in dbid2 (object)
table(gsub(":.*","",crossId$dbid2))

# Generate new table with ids and DBs separaterd in 2 columns
crossId$DB2 <- gsub(":.*","",crossId$dbid2)
crossId$DB1 <- gsub(":.*","",crossId$dbid1)
crossId$id2 <- gsub(".*:","",crossId$dbid2)
crossId$id1 <- gsub(".*:","",crossId$dbid1)
dim(crossId)

## Remove crossIds without a colon (e.g. definitions, ...)
head(grep(":",crossId$dbid1,invert = T,value = T))
head(grep(":",crossId$dbid2,invert = T,value = T))
crossId <- crossId[grepl(":",crossId$dbid2) & grepl(":",crossId$dbid1) ,]
dim(crossId)


## an integer is a correct disease ID > separate those that are (ToKeep) from those that aren't (ToCheck) - warnings are normal
table(!is.na(as.numeric(crossId$id2)))
table(!is.na(as.numeric(crossId$id1)))
toKeep <- crossId[which(!is.na(as.numeric(crossId$id2)) &
                          !is.na(as.numeric(crossId$id1))),]
dim(toKeep)
toCheck <- crossId[-which(!is.na(as.numeric(crossId$id2)) &
                            !is.na(as.numeric(crossId$id1))),]
dim(toCheck)
## When removing prefix (e.g., a letter), an integer is a correct disease ID > change those from ToCheck to ToKeep
table(!is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id2))))
table(!is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id1))))
toKeep <- rbind(toKeep, 
                toCheck[which(!is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id2))) &
                                !is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id1)))),])
dim(toKeep)
toCheck <- toCheck[-which(!is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id2))) &
                            !is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id1)))),]
dim(toCheck)

# Check that DB1 is "EFO" or "Orphanet" (or DOID?)
table(toCheck$DB1)
table(toKeep$DB1)

## In toCheck, check that none of the DB2 are in the list of DBs used for DOD (Data-Source-Model-Repository > 00-Documentation)
# NB there are some DBs with similar names (eg ICD9CM_2006), but these are old versions > don't keep them
table(toCheck$DB2)
toCheck[toCheck$DB2 == "ICD9CM",]
toCheck[toCheck$DB2 == "ICD10CM",]
head(toCheck[toCheck$DB2 == "DO",])

## ## In toKeep
table(toKeep$DB2)
toKeep[grep("url",toKeep$dbid2),]
toKeep[grep("https",toKeep$dbid2),]
toKeep <- toKeep[grep(paste("url","https",sep = "|"),toKeep$dbid2,invert = T),]

# Generate final table (only ToKeep)
crossId <- setNames(toKeep[,c("dbid1","dbid2")],c("id1","id2"))
dim(crossId)
head(crossId)

# Replace some DB names by usual ones
crossId$id1 <- gsub(" ", "", crossId$id1)
crossId$id2 <- gsub(" ", "", crossId$id2)
crossId$id2 <- gsub("MESH","MeSH",crossId$id2)
crossId$id2 <- gsub("\\bNCI\\b","NCIt",crossId$id2)
crossId$id2 <- gsub("ORDO","ORPHA",crossId$id2)
crossId$id2 <- gsub("MEDDRA","MedDRA",crossId$id2)
crossId$id2 <- gsub("UMLS_CUI","UMLS",crossId$id2)
crossId$id2 <- gsub("SNOMECT","SNOMEDCT",crossId$id2)
crossId$id2 <- gsub("SNOMED_CT_US_2018_03_01","SNOMEDCT_US_2018_03_01",crossId$id2)

crossId$id2 <- gsub("ICD-10","ICD10",crossId$id2)
crossId$DB2 <- gsub(":.*","",crossId$id2)
crossId$DB1 <- gsub(":.*","",crossId$id1)
table(crossId$DB2)

## Remove self references
crossId[which(crossId$id1 == crossId$id2),]
## crossId <- crossId[-which(crossId$id1 == crossId$id2),]

######################################
## entryId
# Only select diseases
entryId <- id[id$id %in% disease$descendants,] %>% as_tibble()

# Check all IDs have ":" in them, otherwise remove them
table(gsub(":.*","",entryId$id))
entryId <- entryId[grep(":",entryId$id),,drop = FALSE]
table(gsub(":.*","",entryId$id))
# Generate a DB column, then only select DB, id and def columns 
entryId$DB <- gsub(":.*","",entryId$id)
head(entryId)
entryId <- entryId[,c("DB","id","def")]
# Look for Ids with # in them
unique(grep("#",entryId$id, value =T))
## Empty definition to NA
tail(sort(table(nchar(entryId$def))))
# entryId$def <- ifelse(entryId$def == "",NA,entryId$def)
## Check characters for \t, \n, \r and put to ASCII
#table(unlist(sapply(entryId$def, strsplit, split = "")))
entryId$def <- iconv(x = entryId$def,to="ASCII//TRANSLIT")
entryId$def <- gsub(paste("\n","\t","\r", sep = "|")," ",entryId$def)
entryId$def <- gsub("\"","'",entryId$def)
table(unlist(sapply(entryId$def, strsplit, split = "")))

## Check duplicated records
dim(entryId)
length(unique(entryId[,"id"]))

## check all crossId$id1 are in entryId
table(crossId$id1 %in% entryId$id)

######################################
## idNames

# select diseases only
idNames <- syn[syn$id %in% disease$descendants,] %>%
  as_tibble()
idNames$canonical <- FALSE
# check all ids contain ":"
table(gsub(":.*","",idNames$id))
# check for id with "#" in them
unique(grep("#",idNames$id, value =T))

## Look for labels in id table, check all have ":" and none have "#"
lbl <- id[id$id %in% disease$descendants,c("id","label")] %>% 
  as_tibble()
lbl$canonical <- TRUE
table(gsub(":.*","",lbl$id))
head(lbl)
unique(grep("#",lbl$id, value =T))
# lbl <- lbl[grep("#",lbl$id,invert = T, value = F),]

## 
idNames <- idNames %>%
  as_tibble() %>%
  bind_rows(lbl %>% select(id, syn = label, canonical)) %>%
  mutate(DB = gsub(":.*","", id))
## unique
dim(unique(idNames))
idNames <- idNames[order(idNames$canonical,decreasing = T),]
idNames <- unique(idNames)
dim(idNames)

## Check characters for \t, \n, \r and put to ASCII
# table(unlist(sapply(idNames$syn, strsplit, split = "")))
idNames$syn <- iconv(x = idNames$syn,to="ASCII//TRANSLIT")
idNames$syn <- gsub(paste("\n","\t","\r", sep = "|")," ",idNames$syn)
idNames$syn <- gsub("\"","'",idNames$syn)
idNames$syn <- gsub("\\\\","'",idNames$syn)
table(unlist(sapply(idNames$syn, strsplit, split = "")))

# Remove NA and duplicates
table(is.na(idNames$syn))
idNames <- idNames[!is.na(idNames$syn),]
head(idNames)
idNames <- unique(idNames)
head(idNames)
## check all idNames are in entryId
table(idNames$id %in% entryId$id)

## Remove empty names, ifany + look for problematic names amongst short names
nc <- nchar(idNames$syn)
table(nc)
idNames[which(nc < 3),]
idNames[which(nc == 3),]
idNames[which(nc == 1),]
# idNames <- idNames[-which(nc == 0),]

## Not every ID has a definition available, in this case, the canonical label will be used
entryId <- entryId %>% 
  mutate(def = case_when(is.na(def) ~ lbl$label[match(id,lbl$id)],
                         TRUE ~ def))

######################################
## parentId

# Get edges including diseases, check they all have ":" and generate DB columns
parentId <- edgesJson[which(edgesJson$obj %in% disease$descendants),c("sub","obj")]
table(gsub(":.*","",parentId$sub))
table(gsub(":.*","",parentId$obj))
names(parentId) <- c("id","parent")
parentId$DB <- gsub(":.*","",parentId$id)
parentId$pDB <- gsub(":.*","",parentId$parent)
parentId$origin <- "DO"

## check all parentId are in entryId
table(parentId$id %in% entryId$id)
table(parentId$parent %in% entryId$id)
## "Phenome" not in entryId --> OK
parentId[!(parentId$parent %in% entryId$id),]

## Add levels
getAncestors <- function(id){
  direct <- termParents[[id]]
  parents <- direct
  level <- 0
  dLev <- c()
  for(d in direct){
    dPar <- getAncestors(d)
    dLev <- c(dLev, dPar$level)
    parents <- c(parents, dPar$parents)
  }
  if(length(dLev)>0){
    level <- max(dLev)+1
  }
  return(list(parents=unique(parents), level=level))
}

parentList <- unstack(parentId, parent~id)
termParents <- parentList
library(BiocParallel)
bpparam <- MulticoreParam(workers = 30)

termAncestors <- bplapply(
  parentId$id,
  getAncestors,
  BPPARAM = bpparam
)
names(termAncestors) <- parentId$id

entryId <- entryId %>%
  mutate(
    level=unlist(lapply(termAncestors, function(x) x$level))[entryId$id]
  ) %>%
  mutate(level = case_when(is.na(level) ~ 0,
                           TRUE ~ level))

#######################################
# Alternative ID
altId <- alt[alt$id %in% disease$descendants,] %>%
  as_tibble()
table(altId$id %in% entryId$id)
table(altId$alt %in% entryId$id)
toAdd <- altId %>%
  select(id = alt) %>%
  mutate(DB = "DOID", 
         def = NA, 
         level = NA) %>%
  distinct()
entryId <- entryId %>%
  bind_rows(toAdd)
altId <- altId %>%
  mutate(DB = gsub(":.*", "", id),
         id = gsub(".*:", "", id),
         altDB = gsub(":.*", "", alt),
         alt = gsub(".*:", "", alt)) 

#######################################
# Remove DB from all ID names from all tables
crossId$id1 <- gsub(".*:","",crossId$id1)
crossId$id2 <- gsub(".*:","",crossId$id2)
entryId$id <- gsub(".*:","",entryId$id)
parentId$id <- gsub(".*:","",parentId$id)
parentId$parent <- gsub(".*:","",parentId$parent)
idNames$id <- gsub(".*:","",idNames$id)

############################
# Put order into columns in all tables
DO_idNames <- idNames[,c("DB","id","syn","canonical")]
DO_parentId <- parentId[,c("DB","id","pDB","parent","origin")]
DO_crossId <- crossId[,c("DB1","id1","DB2","id2")]
DO_entryId <- entryId[,c("DB","id","def","level")]
DO_altId <- altId[,c("DB", "id", "altDB", "alt")]

############################
toSave <- grep("^DO[_]", ls(), value=T)
for(f in toSave){
  message(paste("Saving", f))
  ## Ensure unicity
  assign(f, get(f))
  if(length(names(f))==0){
    f <- unique(f)
  }
  ##
  write.table(
    get(f),
    file=file.path(ddir, paste(f, ".txt", sep="")),
    sep="\t",
    row.names=FALSE, col.names=TRUE,
    quote=TRUE,
    qmethod = "double"
  )
}

##############################################################
## Check model
## source("../../00-Utils/autoCheckModel.R")

