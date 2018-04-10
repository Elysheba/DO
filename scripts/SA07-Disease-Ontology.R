rm(list=ls())
gc()

#############################
load("curOptions.rda")
curDate <- Sys.Date()
# curDate <- "2017-05-23"

oboFn <- file.path(tmpDir, paste("doid-", curDate, ".obo", sep=""))
if(!file.exists(oboFn)){
    download.file(
        url="https://github.com/DiseaseOntology/HumanDiseaseOntology/raw/master/src/ontology/doid.obo",
        destfile=oboFn
    )
}
obo <- readLines(oboFn)
obo <- obo[-(which(obo=="[Typedef]")[1]:length(obo))]

############################
## Basic information
starts <- which(obo=="[Term]")
ends <- c(starts[-1]-1, length(obo))
doDef <- do.call(rbind, apply(
    data.frame(starts, ends),
    1,
    function(x){
        termDesc <- obo[(x[1]+1):(x[2]-1)]
        ##
        fn <- "^id: "
        id <- sub(fn, "", grep(fn, termDesc, value=T))
        if(length(id)==0) id <- NA
        ##
        fn <- "^name: "
        name <- sub(fn, "", grep(fn, termDesc, value=T))
        if(length(name)==0) name <- NA
        ##
        fn <- "^def: "
        def <- sub(fn, "", grep(fn, termDesc, value=T))
        def <- sub('".*$', "", sub('^"', "", def))[1]
        if(length(def)==0) def <- NA
        ##
        fn <- "^is_a: "
        parent <- sub(fn, "", grep(fn, termDesc, value=T))
        fn <- " [!].*$"
        parent <- sub(fn, "", parent)
        # if(length(parent)==0) parent <- NA
        parent <- paste(parent, collapse="||")
        ##
        fn <- "^alt_id: "
        altId <- sub(fn, "", grep(fn, termDesc, value=T))
        altId <- paste(unique(c(id, altId)), collapse="||")
        ##
        fn <- "^synonym: "
        synonym <- sub(fn, "", grep(fn, termDesc, value=T))
        synonym <- sub('".*$', "", sub('^"', "", synonym))
        synonym <- paste(unique(c(name, synonym)), collapse="||")
        ##
        fn <- "^xref: EFO:"
        EFO <- sub(fn, "", grep(fn, termDesc, value=T))
        EFO <- paste(EFO, collapse="||")
        ##
        fn <- "^xref: MESH:"
        MESH <- sub(fn, "", grep(fn, termDesc, value=T))
        MESH <- paste(MESH, collapse="||")
        ##
        fn <- "^xref: OMIM:"
        OMIM <- sub(fn, "", grep(fn, termDesc, value=T))
        OMIM <- paste(OMIM, collapse="||")
        ##
        fn <- "^xref: UMLS_CUI:"
        CUI <- sub(fn, "", grep(fn, termDesc, value=T))
        CUI <- paste(CUI, collapse="||")
        ##
        fn <- "^is_obsolete: "
        obsolete <- sub(fn, "", grep(fn, termDesc, value=T))
        if(length(obsolete)==0){
            obsolete <- FALSE
        }else{
            obsolete <- as.logical(obsolete)
        }
        ##
        return(data.frame(
            id=id, name=name, def=def,
            parent=parent,
            altId=altId,
            synonym=synonym,
            EFO=EFO, MESH=MESH, OMIM=OMIM, CUI=CUI,
            obsolete=obsolete,
            stringsAsFactors=F
        ))
    }
))

############################
## Mapping other ID and synonyms
altId <- unique(doDef[, c("id", "altId")])
altIdList <- strsplit(altId$altId, "[|][|]")
names(altIdList) <- altId$id
altId <- stack(altIdList)
colnames(altId) <- c("alt", "id")
altId$id <- as.character(altId$id)
altId$alt <- as.character(altId$alt)
altId <- altId[,c("id", "alt")]
##
synonym <- unique(doDef[, c("id", "synonym")])
synonymList <- strsplit(synonym$synonym, "[|][|]")
names(synonymList) <- synonym$id
synonym <- stack(synonymList)
colnames(synonym) <- c("synonym", "id")
synonym$id <- as.character(synonym$id)
synonym$synonym <- as.character(synonym$synonym)
synonym <- synonym[,c("id", "synonym")]
##
efoMap <- unique(doDef[, c("id", "EFO")])
efoMapList <- strsplit(efoMap$EFO, "[|][|]")
names(efoMapList) <- efoMap$id
efoMap <- stack(efoMapList)
colnames(efoMap) <- c("EFO", "DO")
efoMap$DO <- as.character(efoMap$DO)
efoMap$EFO <- as.character(efoMap$EFO)
efoMap <- efoMap[,c("DO", "EFO")]
##
meshMap <- unique(doDef[, c("id", "MESH")])
meshMapList <- strsplit(meshMap$MESH, "[|][|]")
names(meshMapList) <- meshMap$id
meshMap <- stack(meshMapList)
colnames(meshMap) <- c("MESH", "DO")
meshMap$DO <- as.character(meshMap$DO)
meshMap$MESH <- as.character(meshMap$MESH)
meshMap <- meshMap[,c("DO", "MESH")]
##
omimMap <- unique(doDef[, c("id", "OMIM")])
omimMapList <- strsplit(omimMap$OMIM, "[|][|]")
names(omimMapList) <- omimMap$id
omimMap <- stack(omimMapList)
colnames(omimMap) <- c("OMIM", "DO")
omimMap$DO <- as.character(omimMap$DO)
omimMap$OMIM <- as.character(omimMap$OMIM)
omimMap <- omimMap[,c("DO", "OMIM")]
##
cuiMap <- unique(doDef[, c("id", "CUI")])
cuiMapList <- strsplit(cuiMap$CUI, "[|][|]")
names(cuiMapList) <- cuiMap$id
cuiMap <- stack(cuiMapList)
colnames(cuiMap) <- c("CUI", "DO")
cuiMap$DO <- as.character(cuiMap$DO)
cuiMap$CUI <- as.character(cuiMap$CUI)
cuiMap <- cuiMap[,c("DO", "CUI")]

############################
## Parents
tmp <- unique(doDef[, c("id", "parent")])
doParents <- strsplit(tmp$parent, "[|][|]")
names(doParents) <- tmp$id
doParents <- lapply(
    doParents,
    function(x){
        x <- setdiff(unique(x), NA)
        return(x)
    }
)
getAncestors <- function(doid){
    direct <- doParents[[doid]]
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
doAncestors <- lapply(
    unique(doDef$id),
    getAncestors
)
names(doAncestors) <- names(doParents)
doDef$level <- unlist(lapply(doAncestors, function(x) x$level))[doDef$id]
doAncestors <- lapply(doAncestors, function(x) x$parents)
##
doRel <- rbind(
    stack(doAncestors),
    data.frame(values=names(doAncestors), ind=names(doAncestors))
)
doAncestors <- unstack(doRel, values~ind)
doDescendants <- unstack(doRel, ind~values)

############################
## DO def
doDef <- doDef[, setdiff(
    colnames(doDef),
    c("parent", "altId", "synonym", "EFO", "MESH", "OMIM", "CUI")
)]
rownames(doDef) <- doDef$id

############################
## Disease ontology / HP
fn <- file.path(tmpDir, paste("hp_common_annotations_all-", curDate, ".tab", sep=""))
if(!file.exists(fn)){
    download.file(
        url="http://pubmed-browser.human-phenotype-ontology.org/hp_common_annotations_all.tab",
        destfile=fn
    )
}
doHp <- read.table(
    file=fn,
    header=FALSE,
    quote="", comment="",
    sep="\t",
    stringsAsFactors=F
)[,c(1, 2 , 3, 4)]
colnames(doHp) <- c("MESH", "Name", "DO", "HP")

hpMesh <- data.frame(
    "MESH"=doHp$MESH,
    "name"=doHp$Name,
    HP=doHp$HP,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

meshNotInDo <- setdiff(doHp$MESH, meshMap$MESH)
meshNotInDo <- unique(doHp[
    which(doHp$MESH %in% meshNotInDo),
    c("MESH", "Name", "DO")
    ])
meshNotInDo$sdo <- sub("DOID:0*", "", meshNotInDo$DO)
sdoToDo <- data.frame(
    RealDO=doDef$id,
    DOName=doDef$name,
    obsolete=doDef$obsolete,
    stringsAsFactors=FALSE
)
rownames(sdoToDo) <- sub("DOID:0*", "", doDef$id)
meshNotInDo <- cbind(meshNotInDo, sdoToDo[meshNotInDo$sdo,])
toAdd <- meshNotInDo[which(meshNotInDo$DO!=""),c("RealDO", "MESH")]
colnames(toAdd) <- c("DO", "MESH")
meshMap <- rbind(meshMap, toAdd)

hpDo <- doHp[which(doHp$MESH %in% meshMap$MESH), c("MESH", "HP")]
hpDo <- merge(
    hpDo,
    meshMap,
    by="MESH",
    all.x=T,
    all.y=F
)
hpDo <- unique(data.frame(
    DO=hpDo$DO,
    HP=hpDo$HP,
    stringsAsFactors=FALSE,
    check.names=FALSE
))

############################
## Save
do.definitions <- doDef
do.parents <- doParents
do.ancestors <- doAncestors
do.descendants <- doDescendants 
do.altid <- altId
do.synonym <- synonym
##
do.efo <- efoMap
do.mesh <- meshMap
do.omim <- omimMap
do.cui <- cuiMap
##
do.hp <- hpDo
mesh.hp <- hpMesh
##
toSave <- c(
    grep("^do[.]", ls(), value=T),
    grep("^mesh[.]", ls(), value=T)
)
for(f in toSave){
    message(paste("Saving", f))
    save(
        list=f,
        file=file.path(tmpDir, paste(f, ".rda", sep=""))
    )
}
