## FIXME: we should offer to determine packageName based on header of
## CDF file and oligo::cleanPlatformName.

setMethod("makePdInfoPackage", "AffySNPPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
              validInput <- function(x, destPath) {
                  msg <- NULL
                  ok <- sapply(c("cdfFile", "csvAnnoFile", "csvSeqFile"),
                               function(slt) file.exists(slot(x, slt)))
                  if (!all(ok))
                    msg <-
                      paste("missing file(s):",
                            paste(sapply(names(ok[!ok]), function(slt) slt),
                                  "='",
                                  sapply(names(ok[!ok]),
                                         function(slt) slot(x, slt)),
                                  "'",
                                  collapse=", ", sep=""))
                  if (file.exists(destPath))
                    msg <-
                      c(msg,
                        paste("destination exists, remove or rename: '",
                              destPath, "'", sep=""))
                  if (is.null(msg)) TRUE else msg
              }
              chip <- chipName(object)
              pkgName <- cleanPlatformName(chip)
              valid <- validInput(object, file.path(destDir, pkgName))
              if (!identical(valid, TRUE))
                stop(paste(valid, collapse="\n  "))
              extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
              dbFileName <- paste(pkgName, "sqlite", sep=".")
              dbFilePath <- file.path(extdataDir, dbFileName)
              seqMatFile <- file.path(extdataDir, "seqMat.rda")
              geometry <- paste(readCdfHeader(object@cdfFile)[c("nrows", "ncols")], collapse=";")
              syms <- list(MANUF=object@manufacturer,
                           VERSION=object@version,
                           GENOMEBUILD=object@genomebuild,
                           AUTHOR=object@author,
                           AUTHOREMAIL=object@email,
                           LIC=object@license,
                           DBFILE=dbFileName,
                           CHIPNAME=chip,
                           PKGNAME=pkgName,
                           PDINFONAME=pkgName,
                           PDINFOCLASS="AffySNPPDInfo",
                           GEOMETRY=geometry)

              templateDir <- system.file("pd.PKG.template",
                                         package="pdInfoBuilder")
              createPackage(pkgname=pkgName, destinationDir=destDir,
                            originDir=templateDir, symbolValues=syms,
# getting error on invalid use of ..., vjc 5 jul 2007
#                            quiet=quiet, ...)
                            quiet=quiet)
              dir.create(extdataDir, recursive=TRUE)
              buildPdInfoDb(object@cdfFile, object@csvAnnoFile,
                            object@csvSeqFile, dbFilePath, seqMatFile,
                            batch_size=batch_size, verbose=!quiet)
              ## copy external resource files
              extFiles <- c(object@splineParamFile,
                            object@crlmmInfoFile,
                            object@referenceDistFile)
              for (ef in extFiles) {
                  if (nchar(ef) > 0)
                    file.copy(ef, extdataDir)
              }
          })

setMethod("makePdInfoPackage", "AffySNPCNVPDInfoPkgSeed",
          function(object, destDir=".", batch_size=1000, quiet=FALSE, unlink=FALSE) {
              chip <- chipName(object)
              pkgName <- cleanPlatformName(chip)
              extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
              dbFileName <- paste(pkgName, "sqlite", sep=".")
              dbFilePath <- file.path(extdataDir, dbFileName)
              seqMatFile <- file.path(extdataDir, "seqMat.rda")
              seqMatFileCnv <- file.path(extdataDir, "seqMatCNV.rda")
              geometry <- paste(readCdfHeader(object@cdfFile)[c("nrows", "ncols")], collapse=";")
              syms <- list(MANUF=object@manufacturer,
                           VERSION=object@version,
                           GENOMEBUILD=object@genomebuild,
                           AUTHOR=object@author,
                           AUTHOREMAIL=object@email,
                           LIC=object@license,
                           DBFILE=dbFileName,
                           CHIPNAME=chip,
                           PKGNAME=pkgName,
                           PDINFONAME=pkgName,
                           PDINFOCLASS="AffySNPCNVPDInfo",
                           GEOMETRY=geometry)

              templateDir <- system.file("pd.PKG.template",
                                         package="pdInfoBuilder")
              createPackage(pkgname=pkgName, destinationDir=destDir,
                            originDir=templateDir, symbolValues=syms,
                            quiet=quiet)
              dir.create(extdataDir, recursive=TRUE)
              snp6.buildPdInfoDb(object@cdfFile, object@csvAnnoFile,
                                 object@csvSeqFile, object@csvAnnoFileCnv,
                                 object@csvSeqFileCnv, dbFilePath, seqMatFile, seqMatFileCnv,
                                 batch_size=batch_size, verbose=!quiet)
              ## copy external resource files
              extFiles <- c(object@splineParamFile,
                            object@crlmmInfoFile,
                            object@referenceDistFile)
              for (ef in extFiles) {
                  if (nchar(ef) > 0)
                    file.copy(ef, extdataDir)
              }
          })
          
## modified by Matt Settles June 2,2008
setMethod("makePdInfoPackage", "NgsExpressionPDInfoPkgSeed",
    function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
  
## debug
#object = pkg
#destDir=targ
#batch_size=10000
#quiet=FALSE
## end debug
    validInput <- function(x, destPath) {
        msg <- NULL
        ok <- sapply(c("pairFile", "xysFile"),
            function(slt) file.exists(slot(x, slt)))
        if (all(!ok))
            msg <- "Must specify one of: pairFile or xysFile"
        else if(all(ok))
            msg <- "Must specify ONLY one of: pairFile or xysFile"
        ok <- sapply(c("ndfFile"),
            function(slt) file.exists(slot(x, slt)))
        if (!all(ok))
            msg <- paste(msg,"\n","missing file(s):",
                paste(sapply(names(ok[!ok]), function(slt) slt),
                    "='",
                    sapply(names(ok[!ok]), function(slt) slot(x, slt)),
                    "'",
                    collapse=", ", sep=""))
        if (file.exists(destPath))
            msg <- c(msg,"\n",
                paste("destination exists, remove or rename: '",
                    destPath, "'", sep=""))
        if (is.null(msg)) TRUE else msg
    }

    chip <- chipName(object) ## only place that uses xys or pair file
    pkgName <- cleanPlatformName(chip)
    valid <- validInput(object, file.path(destDir, pkgName))
    if (!identical(valid, TRUE))
        stop(paste(valid, collapse="\n  "))
    extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
    dbFileName <- paste(pkgName, "sqlite", sep=".")
    dbFilePath <- file.path(extdataDir, dbFileName)
    ## get geometry from ndf File - MS
    ndfdata <- read.delim(object@ndfFile, as.is=TRUE, header=TRUE)
    geometry <- paste(max(ndfdata$Y), max(ndfdata$X), sep=";")
    nx <- max(ndfdata$X)
    
    ## I don't see the point to all below
    #ndf.idx <- as.integer(ndfdata$X+(ndfdata$Y-1)*nx)
    ### do we really need any of this
    #if(file.exists(object@pairFile)) { # pairFile used
    #    pairdata <- read.delim(object@pairFile, as.is=TRUE, header=TRUE, comment.char="#")
    #    pair.idx <- as.integer(pairdata$X+(pairdata$Y-1)*nx)
    #    idx <- match(pair.idx, ndf.idx) ## match ndf to pair, removes Control probes 
    #    ndfdata <- ndfdata[idx,]
    #    rm(pairdata,pair.idx,idx)
    #} else if (file.exists(object@xysfile)) { #xysFile used
    #    xysdata <- read.delim(object@xysFile, as.is=TRUE, header=TRUE, comment.char="#")
    #    xys.idx <- as.integer(xysdata$X+(xysdata$Y-1)*nx)
    #    idx <- match(xys.idx, ndf.idx) ## match ndf to pair              
    #    ndfdata <- ndfdata[idx,]
    #    idx <- which(is.na(xysdata$SIGNAL)) 
    #    if(length(idx) > 0) ndfdata[idx, "CONTAINER"] <- "OLIGO_CONTROL" # added if statement - MS                
    #    rm(xysdata,xys.idx,idx) 
    #}
    ## in case one removes  H_CODE,  NGS_CONTROLS and V_CODE from the NDF file
    ## in case two renames  H_CODE,  NGS_CONTROLS and V_CODE to OLIGO_CONTROL
    # rm(ndf.idx, idx);gc();gc();

    ndfdata <- cbind(fid=as.integer(ndfdata$X+(ndfdata$Y-1)*nx) , ndfdata) 
    ndfdata <- as.data.frame(ndfdata)
    ndfdata <- ndfdata[order(ndfdata[["fid"]]),]
    ndfdata[["fid"]] <- 1:nrow(ndfdata)  
    syms <- list(MANUF=object@manufacturer,
        VERSION=object@version,
        GENOMEBUILD=object@genomebuild,
        AUTHOR=object@author,
        AUTHOREMAIL=object@email,
        LIC=object@license,
        DBFILE=dbFileName,
        CHIPNAME=chip,
        PKGNAME=pkgName,
        PDINFONAME=pkgName,
        PDINFOCLASS="NgsExpressionPDInfo",
        GEOMETRY=geometry)
    templateDir <- system.file("pd.PKG.template", package="pdInfoBuilder")
    createPackage(pkgname=pkgName, destinationDir=destDir,
    originDir=templateDir, symbolValues=syms, quiet=quiet)
    dir.create(extdataDir, recursive=TRUE)
    
    buildPdInfoDb.ngsExprs(ndfdata,object@ngdFile, dbFilePath, # changed to buildPdInfoDb.ngsExprs - MS
        batch_size=batch_size, verbose=!quiet)
})

## TODO: fix for tiling arrays
## modified by Matt Settles June 2,2008
setMethod("makePdInfoPackage", "NgsTilingPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
            ## the validInput that was once here was moved to setValidity
            chip <- chipName(object)
            pkgName <- cleanPlatformName(chip)
            extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
            dbFileName <- paste(pkgName, "sqlite", sep=".")
            dbFilePath <- file.path(extdataDir, dbFileName)
            ## without geometry in header file this is the best alternative - MS
            ## FIXME!!! need geometry 
            ## Awww... it would be great if NGS had a header describing the array.
            ## ndfdata was supposed to appear only on loadByBatch
            ## appears here to avoid reading the same file twice
            ndfdata <- read.delim(object@ndfFile, as.is=TRUE, header=TRUE)
            geometry <- paste(max(ndfdata$Y), max(ndfdata$X), sep=";")
            nx <- max(ndfdata$X)
            ndf.idx <- as.integer(ndfdata$X+(ndfdata$Y-1)*nx)
            if(file.exists(object@pairFile)) { # pairFile used
              pairdata <- read.delim(object@pairFile, as.is=TRUE, header=TRUE, comment.char="#")
              pair.idx <- as.integer(pairdata$X+(pairdata$Y-1)*nx)
              idx <- match(pair.idx, ndf.idx) ## match ndf to pair
              ndfdata <- ndfdata[idx,]
              idx <- which(is.na(pairdata$PM))
              rm(pairdata,pair.idx) 
            }else if (file.exists(object@xysFile)) { #xysFile used
              xysdata <- read.delim(object@xysFile, as.is=TRUE, header=TRUE, comment.char="#")
              xys.idx <- as.integer(xysdata$X+(xysdata$Y-1)*nx)
              idx <- match(xys.idx, ndf.idx) ## match ndf to xys
              ndfdata <- ndfdata[idx,]
              idx <- which(is.na(xysdata$SIGNAL))
              rm(xysdata,xys.idx) 
            }
            if(length(idx) > 0) ndfdata[idx, "CONTAINER"] <- "OLIGO_CONTROL" # added if statement - MS                
            ndfdata <- cbind(fid=as.integer(ndfdata$X+(ndfdata$Y-1)*nx), ndfdata)
            ndfdata <- as.data.frame(ndfdata)
            ndfdata <- ndfdata[order(ndfdata[["fid"]]),]
            ndfdata[["fid"]] <- 1:nrow(ndfdata)
            rm(ndf.idx, idx);gc();gc()
            syms <- list(MANUF=object@manufacturer,
                         VERSION=object@version,
                         GENOMEBUILD=object@genomebuild,
                         AUTHOR=object@author,
                         AUTHOREMAIL=object@email,
                         LIC=object@license,
                         DBFILE=dbFileName,
                         CHIPNAME=chip,
                         PKGNAME=pkgName,
                         PDINFONAME=pkgName,
                         PDINFOCLASS="NgsTilingPDInfo",
                         GEOMETRY=geometry)
            templateDir <- system.file("pd.PKG.template",
                                       package="pdInfoBuilder")
            createPackage(pkgname=pkgName, destinationDir=destDir,
                          originDir=templateDir, symbolValues=syms,
                          quiet=quiet)
            dir.create(extdataDir, recursive=TRUE)
            buildPdInfoDb.ngsTiling(ndfdata, object@posFile, dbFilePath, ## changed to buildPdInfoDb.ngsTiling - MS
                                    batch_size=batch_size, verbose=!quiet) 
          })


setMethod("makePdInfoPackage", "AffyExpressionPDInfoPkgSeed",
    function(object, destDir=".", batch_size=10000, quiet=FALSE, unlink=FALSE) {
        valid <- validInput(object,c("cdfFile"), c("csvAnnoFile", "tabSeqFile"))
        if (!identical(valid, TRUE))
            stop(paste(valid, collapse="\n  "))

        if(is.na(object@chipName)) object@chipName <- chipName(object)
        pkgName <- cleanPlatformName(object@chipName)
        
        dbFileName <- paste(pkgName, "sqlite", sep=".")        
        setupPackage(object,pkgName,destDir,dbFileName,unlink,quiet)
        
        extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
        dbFilePath <- file.path(extdataDir, dbFileName)
        seqMatFile <- file.path(extdataDir, "seqMat.rda")
        
        buildPdInfoDb.affyExpr(
                object@cdfFile, 
                object@csvAnnoFile,
                object@tabSeqFile, 
                dbFilePath, 
                seqMatFile,
                batch_size=batch_size, 
                verbose=!quiet)
})

          
setMethod("makePdInfoPackage", "AffyTilingPDInfoPkgSeed",
    function(object, destDir=".", batch_size=1, quiet=FALSE, unlink=FALSE) {
## debug
#destDir="."; batch_size=1; quiet=FALSE; unlink=TRUE
## end debug
        valid <- validInput(object,c("bpmapFile", "celFile"))
        if (!identical(valid, TRUE))
            stop(paste(valid, collapse="\n  "))
              
        if(is.na(object@chipName)) object@chipName <- chipName(object)
        pkgName <- cleanPlatformName(object@chipName)
              
        dbFileName <- paste(pkgName, "sqlite", sep=".")        
        setupPackage(object,pkgName,destDir,dbFileName,unlink,quiet)
        
        extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
        dbFilePath <- file.path(extdataDir, dbFileName)
        seqMatFile <- file.path(extdataDir, "seqMat.rda")
              
        buildPdInfoDb.affyTiling(
                object@bpmapFile, 
                object@celFile,
                dbFilePath, 
                seqMatFile,
                batch_size=batch_size, 
                verbose=!quiet)
})


setMethod("makePdInfoPackage", "AffySTPDInfoPkgSeed",
    function(object, destDir=".", batch_size=1, quiet=FALSE, unlink=FALSE) {
## debug
#destDir="."; batch_size=1; quiet=FALSE; unlink=TRUE
## end debug

        valid <- validInput(object,c("pgfFile", "clfFile"),c("probeFile", "transFile"))
        if (!identical(valid, TRUE))
            stop(paste(valid, collapse="\n  "))
      
        if(is.na(object@chipName)) object@chipName <- chipName(object)
        pkgName <- cleanPlatformName(object@chipName)
      
        dbFileName <- paste(pkgName, "sqlite", sep=".")        
        setupPackage(object,pkgName,destDir,dbFileName,unlink,quiet)
            
        extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
        dbFilePath <- file.path(extdataDir, dbFileName)
        seqMatFile <- file.path(extdataDir, "seqMat.rda")
        
#        buildPdInfoDb.affyST(
#                object@pgfFile, 
#                object@clfFile,
#                object@probeFile,
#                object@transFile, 
#                dbFilePath,
#                seqMatFile,
#                batch_size=batch_size,
#                verbose=!quiet)
            clf <- readClfEnv(object@clfFile)
            pgf <- readPgfEnv(object@pgfFile)
            
            db <- initDb.affyST(dbFilePath)
            t <- ST(loadUnits.affyST(db, pgf, clf))
            if (!quiet) printTime("loadUnitsByBatch", t[3])
            if(!is.na(slot(object,"transFile"))){
                t <- ST(loadAffyCsv.affyST(db, object@transFile, batch_size=batch_size))
                if (!quiet) printTime("loadAffyCsv", t[3])
            }
            if(!is.na(slot(object,"probeFile"))){
            t <- ST(loadAffySeqCsv.affyST(db, object@probeFile, pgf, batch_size=batch_size))
            if (!quiet) printTime("loadAffySeqCsv", t[3])
            }
            t <- ST({
                        sortFeatureTables.affyST(db)
                        createIndicesDb.affyST(db)
                        createTableInfoTable.affyST(db)
                        createFeatureTableInfo.affyST(db)
                    })
            if (!quiet) printTime("DB sort, index creation", t[3])
            
            ##     t <- ST({
            ##         seqMat <- createSeqMat(db)
            ##         save(seqMat, file=matFile, compress=TRUE)
            ##     })
            ##     printTime("sequence matrix", t[3])
            closeDb(db)
})

