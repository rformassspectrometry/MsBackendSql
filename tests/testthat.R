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

mm14_file <- system.file("microtofq", "MM14.mzML", package = "msdata")
mm_db <- dbConnect(SQLite(), tempfile())
createMsqlBackendDatabase(mm_db, c(mm8_file, mm14_file))
mm_be <- backendInitialize(MsqlBackend(), mm_db)

tmt_file <- proteomics(full.names = TRUE)[4L]
tmt_mzr <- backendInitialize(MsBackendMzR(), tmt_file)
tmt_db <- dbConnect(SQLite(), tempfile())
createMsqlBackendDatabase(tmt_db, tmt_file)
tmt_be <- backendInitialize(MsqlBackend(), tmt_db)

test_check("MsqlBackend")

dbDisconnect(mm8_db)
