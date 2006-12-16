setClass("PDInfoPkgSeed",
         representation=representation(
           chipName="character",
           manufacturer="character",
           pdInfoObjectName="character"))

setClass("AffySNPPDInfoPkgSeed",
         contains="PDInfoPkgSeed",
         representation=representation(
           cdfFile="character",
           csvAnnoFile="character",
           csvSeqFile="character"
           ),
         prototype=list(manufacturer="Affymetrix"))

##pkg <- new("AffySNPPDInfoPkgSeed", chipName="Mapping250K_Nsp", pdInfoObjectName="mapping250k.nsp", cdfFile="", csvAnnoFile="", csvSeqFile="")
##makeDataPackage(pkg, "Seth Falcon", "sfalcon@fhcrc.org", "pd.Mapping250K.Nsp", "0.0.9", biocViews="AnnotationData", filePath=".")

setMethod(Biobase::makeDataPackage,
          signature(object="AffySNPPDInfoPkgSeed"),
          function(object, author, email, packageName, packageVersion, license, biocViews,
                   filePath, quiet=FALSE) {

              if (missing(email) || !is.character(email)
                  || length(email) != 1 || grep("@", email) != 1)
                stop("invalid email address")

              extdataDir <- file.path(filePath, packageName, "inst", "extdata")
              dbFileName <- paste(packageName, "sqlite", sep=".")
              dbFilePath <- file.path(extdataDir, dbFileName)
              seqMatFile <- file.path(extdataDir, "seqMat.rda")
              syms <- list(MANUF=object@manufacturer,
                           VERSION=packageVersion,
                           AUTHOR=author,
                           AUTHOREMAIL=email,
                           LIC=license,
                           DBFILE=dbFileName,
                           CHIPNAME=object@chipName,
                           PKGNAME=packageName,
                           PDINFONAME=object@pdInfoObjectName,
                           PDINFOCLASS="AffySNPPDInfo")

              templateDir <- system.file("pd.PKG.template",
                                         package="pdInfoBuilder")
              createPackage(pkgname=packageName, destinationDir=filePath,
                            originDir=templateDir,
                            symbolValues=syms, quiet=quiet)
              dir.create(extdataDir, recursive=TRUE)
              buildPdInfoDb(object@cdfFile, object@csvAnnoFile,
                            object@csvSeqFile, dbFilePath, seqMatFile)
          })
