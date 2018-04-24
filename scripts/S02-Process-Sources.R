setwd("~/Shared/Data-Science/Data-Source-Model-Repository/DO/scripts/")

library(XML)
library(parallel)
library(git2r)
library(RJSONIO)

##
mc.cores <- 55
sdir <- "../sources/DO/src/ontology/"
ddir <- "../data"

###############################################################################@
## Source information ----
###############################################################################@

desc <- RJSONIO::readJSONStream("../DESCRIPTION.json")

sourceFiles <- desc$"source files"
sfi_name <- unlist(lapply(
  sourceFiles,
  function(sf){
    toRet <- sf$"name"
    return(toRet)
  }
))

###############################################################################@
## Data from doid.json ----
###############################################################################@
readJson <- jsonlite::fromJSON(txt = file.path(sdir,"doid.json"))

###########################################
## nodes (id, def, name, xref, label)
nodesJson <- do.call(rbind,
                     lapply(1:nrow(readJson$graphs$nodes[[1]]),
                            function(i){
                              ## id
                              id <- gsub("_",":",gsub(".*/","",readJson$graphs$nodes[[1]]$id[[i]]))
                              ## def
                              def <- readJson$graphs$nodes[[1]]$meta$definition[i, "val"]
                              ## Xref
                              Xref <- paste(readJson$graphs$nodes[[1]]$meta$xrefs[[i]]$val,collapse = ", ")
                              ## name
                              name <- paste(readJson$graphs$nodes[[1]]$meta$synonyms[[i]]$val,collapse = ", ")
                              ## lbl
                              lbl <- readJson$graphs$nodes[[1]]$lbl[[i]]
                              ## df
                              df <- data.frame(
                                id = id,
                                Xref = Xref,
                                def = def,
                                name = name,
                                label = lbl,
                                stringsAsFactors = FALSE)
                              return(df)
                            }
                     )
)

## edges (parents)
edgesJson <- readJson$graphs$edges[[1]]
edgesJson <- edgesJson[which(edgesJson$pred %in% c("is_a")),]
edgesJson$sub <- gsub("_",":",gsub(".*/","",edgesJson$sub))
edgesJson$obj <- gsub("_",":",gsub(".*/","",edgesJson$obj))

######################################
## crossId
crossId <- unique(nodesJson[,c("id","Xref")])
crossIdList <- strsplit(crossId$Xref, split = ",")
names(crossIdList) <- crossId$id
crossId <- stack(crossIdList)
names(crossId) <- c("id2","id1")
crossId$id2 <- gsub("MESH","MeSH",crossId$id2)
crossId$id2 <- gsub("NCiT","NCIt",crossId$id2)
crossId$id2 <- gsub("NCI","NCIt",crossId$id2)
crossId$id2 <- gsub("ORDO","ORPHA",crossId$id2)
crossId$id2 <- gsub("UMLS_CUI","MedGen",crossId$id2)
crossId$id2 <- gsub("UMLS","MedGen",crossId$id2)
crossId$DB2 <- gsub(":.*","",crossId$id2)
crossId$DB1 <- gsub(":.*","",crossId$id1)

######################################
## entryId
entryId <- nodesJson["id"]
entryId$DB <- gsub(":.*","",entryId$id)
entryId <- entryId[,c("DB","id")]
entryId$definition <- nodesJson$def[match(entryId$id,nodesJson$id)]

######################################
## idNames
idNames <- unique(nodesJson[,c("id","name")])
idNamesList <- strsplit(idNames$name, split = ",")
names(idNamesList) <- idNames$id
idNames <- stack(idNamesList)
names(idNames) <- c("name","id")
## Labels
lbl <- unique(nodesJson[,c("label","id")])
## 
idNames <- rbind(idNames,setNames(lbl, nm = names(idNames)))
idNames$DB <- gsub(":.*","",idNames$id)
idNames$canonical <- ifelse(idNames$name %in% lbl$label, TRUE, FALSE)

######################################
## parentId
parentId <- edgesJson[,c("sub","obj")]
names(parentId) <- c("id","parent")
parentId$DB <- gsub(":.*","",parentId$id)
parentId$pDB <- gsub(":.*","",parentId$parent)
parentId <- parentId[parentId$id %in% nodesJson$id,]
parentId <- parentId[parentId$parent %in% nodesJson$id,]

#######################################
crossId$id1 <- gsub(".*:","",crossId$id1)
crossId$id2 <- gsub(".*:","",crossId$id2)
entryId$id <- gsub(".*:","",entryId$id)
parentId$id <- gsub(".*:","",parentId$id)
parentId$parent <- gsub(".*:","",parentId$parent)
idNames$id <- gsub(".*:","",idNames$id)

############################
DO_idNames <- idNames[,c("DB","id","name","canonical")]
DO_parentId <- parentId[,c("DB","id","pDB","parent")]
DO_crossId <- crossId[,c("DB1","id1","DB2","id2")]
DO_entryId <- entryId[,c("DB","id","definition")]

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
    quote=FALSE
  )
}


