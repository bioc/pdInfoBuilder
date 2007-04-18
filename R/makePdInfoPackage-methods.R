## FIXME: we should offer to determine packageName based on header of
## CDF file and oligo::cleanPlatformName.

setMethod("makePdInfoPackage", "AffySNPPDInfoPkgSeed",
          function(object, destDir=".", batch_size=10000, quiet=FALSE) {
              chip <- chipName(object)
              pkgName <- cleanPlatformName(chip)
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
                            quiet=quiet)
              dir.create(extdataDir, recursive=TRUE)
              buildPdInfoDb(object@cdfFile, object@csvAnnoFile,
                            object@csvSeqFile, dbFilePath, seqMatFile,
                            batch_size=batch_size, verbose=!quiet)
          })
