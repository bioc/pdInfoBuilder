## FIXME: we should offer to determine packageName based on header of
## CDF file and oligo::cleanPlatformName.

setMethod("makePdInfoPackage", "AffySNPPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE) {
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
          function(object, destDir=".", batch_size=10000, quiet=FALSE) {
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

setMethod("makePdInfoPackage", "NgsExpressionPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE) {
              chip <- chipName(object)
              pkgName <- cleanPlatformName(chip)
              extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
              dbFileName <- paste(pkgName, "sqlite", sep=".")
              dbFilePath <- file.path(extdataDir, dbFileName)
              ## FIXME!!! need geometry
              ## Awww... it would be great if NGS had a header describing the array.
              ## ndfdata was supposed to appear only on loadByBatch
              ## appears here to avoid reading the same file twice
              ndfdata <- read.delim(object@ndfFile, as.is=TRUE, header=TRUE)
              xysdata <- read.delim(object@xysFile, as.is=TRUE, header=TRUE, comment.char="#")
              geometry <- paste(max(ndfdata$Y), max(ndfdata$X), sep=";")
              nx <- max(ndfdata$X)
              ndf.idx <- as.integer(ndfdata$X+(ndfdata$Y-1)*nx)
              xys.idx <- as.integer(xysdata$X+(xysdata$Y-1)*nx)
              idx <- match(xys.idx, ndf.idx)
              ndfdata <- ndfdata[idx,]
              idx <- which(is.na(xysdata$SIGNAL))
              ndfdata[idx, "CONTAINER"] <- "OLIGO_CONTROL"
              ndfdata <- cbind(fid=as.integer(ndfdata$X+(ndfdata$Y-1)*nx), ndfdata)
              ndfdata <- as.data.frame(ndfdata)
              ndfdata <- ndfdata[order(ndfdata[["fid"]]),]
              ndfdata[["fid"]] <- 1:nrow(ndfdata)
              rm(xysdata, ndf.idx, xys.idx, idx)
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

              templateDir <- system.file("pd.PKG.template",
                                         package="pdInfoBuilder")
              createPackage(pkgname=pkgName, destinationDir=destDir,
                            originDir=templateDir, symbolValues=syms,
                            quiet=quiet)
              dir.create(extdataDir, recursive=TRUE)
              buildPdInfoDb.ngs(ndfdata, object@posFile, dbFilePath,
                                batch_size=batch_size, verbose=!quiet)
          })

setMethod("makePdInfoPackage", "NgsTilingPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE) {
              chip <- chipName(object)
              pkgName <- cleanPlatformName(chip)
              extdataDir <- file.path(destDir, pkgName, "inst", "extdata")
              dbFileName <- paste(pkgName, "sqlite", sep=".")
              dbFilePath <- file.path(extdataDir, dbFileName)
              ## FIXME!!! need geometry
              ## Awww... it would be great if NGS had a header describing the array.
              ## ndfdata was supposed to appear only on loadByBatch
              ## appears here to avoid reading the same file twice
              ndfdata <- read.delim(object@ndfFile, as.is=TRUE, header=TRUE)
              xysdata <- read.delim(object@xysFile, as.is=TRUE, header=TRUE, comment.char="#")
              geometry <- paste(max(ndfdata$Y), max(ndfdata$X), sep=";")
              nx <- max(ndfdata$X)
              ndf.idx <- as.integer(ndfdata$X+(ndfdata$Y-1)*nx)
              xys.idx <- as.integer(xysdata$X+(xysdata$Y-1)*nx)
              idx <- match(xys.idx, ndf.idx)
              ndfdata <- ndfdata[idx,]
              idx <- which(is.na(xysdata$SIGNAL))
              ndfdata[idx, "CONTAINER"] <- "OLIGO_CONTROL"
              ndfdata <- cbind(fid=as.integer(ndfdata$X+(ndfdata$Y-1)*nx), ndfdata)
              ndfdata <- as.data.frame(ndfdata)
              ndfdata <- ndfdata[order(ndfdata[["fid"]]),]
              ndfdata[["fid"]] <- 1:nrow(ndfdata)
              rm(xysdata, ndf.idx, xys.idx, idx)
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
              buildPdInfoDb.ngs(ndfdata, object@posFile, dbFilePath,
                                batch_size=batch_size, verbose=!quiet)
          })


setMethod("makePdInfoPackage", "AffyExpressionPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE) {
              validInput <- function(x, destPath) {
                  msg <- NULL
                  ok <- sapply(c("cdfFile", "csvAnnoFile", "tabSeqFile"),
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
                           PDINFOCLASS="AffyExpressionPDInfo",
                           GEOMETRY=geometry)

              templateDir <- system.file("pd.PKG.template",
                                         package="pdInfoBuilder")
              createPackage(pkgname=pkgName, destinationDir=destDir,
                            originDir=templateDir, symbolValues=syms,
                            quiet=quiet)
              dir.create(extdataDir, recursive=TRUE)
              buildPdInfoDb.affyExpr(object@cdfFile, object@csvAnnoFile,
                                     object@tabSeqFile, dbFilePath, seqMatFile,
                                     batch_size=batch_size, verbose=!quiet)
          })
