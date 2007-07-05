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
