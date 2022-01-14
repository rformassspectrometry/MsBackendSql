library(testthat)
library(MsqlBackend)
library(Spectra)
library(RSQLite)
library(msdata)

mm8_file <- system.file("microtofq", "MM8.mzML", package = "msdata")
mm8_sps <- Spectra(mm8_file)
mm8_db <- dbConnect(SQLite(), tempfile())
createMsqlBackendDatabase(mm8_db, mm8_file)
mm8_be <- backendInitialize(MsqlBackend(), mm8_db)

test_check("MsqlBackend")

dbDisconnect(mm8_db)
