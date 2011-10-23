library("@PKGNAME@")
library("RUnit")

options(warn=1)

## RUnit Test Suites

dirs <- '.'
testFilePat <- ".*_test\\.R$"

allSuite <- defineTestSuite(name="Test Suite for @PKGNAME@",
                            dirs=dirs,
                            testFileRegexp=testFilePat,
                            rngKind="default",
                            rngNormalKind="default")

testData <- runTestSuite(allSuite)
printTextProtocol(testData, showDetails=FALSE)

