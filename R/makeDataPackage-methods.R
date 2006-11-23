setClass("PDInfoPkg",
         representation=representation(
           chipName="character",
           manufacturer="character",
           pdInfoObjectName="character"))

setClass("PDInfoAffySNPPkg",
         contains="PDInfoPkg",
         representation=representation(
           cdfFile="character",
           csvAnnoFile="character",
           csvSeqFile="character"
           ))

##pkg <- new("PDInfoAffySNPPkg", chipName="Mapping250nsp", manufacturer="Affymetrix", pdInfoObjectName="map250nsp", cdfFile="", csvAnnoFile="", csvSeqFile="")
##makeDataPackage(pkg, "Seth Falcon", "sfalcon@fhcrc.org", "pd.Map250nsp", "0.0.7", biocViews="AnnotationData", filePath=".")

setMethod(Biobase::makeDataPackage,
          signature(object="PDInfoAffySNPPkg"),
          function(object, author, email, packageName, packageVersion, license, biocViews,
                   filePath, quiet=FALSE) {

              if (missing(email) || !is.character(email)
                  || length(email) != 1 || grep("@", email) != 1)
                stop("invalid email address")

              extdataDir <- file.path(filePath, packageName, "inst", "extdata")
              dbFileName <- paste(packageName, "sqlite", sep=".")
              dbFilePath <- file.path(extdataDir, dbFileName)
              syms <- list(MANUF=object@manufacturer,
                           VERSION=packageVersion,
                           AUTHOR=author,
                           AUTHOREMAIL=email,
                           LIC=license,
                           DBFILE=dbFileName,
                           CHIPNAME=object@chipName,
                           PKGNAME=packageName,
                           PDINFONAME=object@pdInfoObjectName,
                           PDINFOCLASS="SNPPDInfo")

              templateDir <- system.file("pd.PKG.template",
                                         package="pdInfoBuilder")
              createPackage(pkgname=packageName, destinationDir=filePath,
                            originDir=templateDir,
                            symbolValues=syms, quiet=quiet)
              dir.create(extdataDir, recursive=TRUE)
              buildPdInfoDb(object@cdfFile, object@csvAnnoFile,
                            object@csvSeqFile, dbFilePath)
          })
