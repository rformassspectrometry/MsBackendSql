library(testthat)
library(MsBackendSql)
library(Spectra)
library(RSQLite)
library(msdata)

mm8_file <- system.file("microtofq", "MM8.mzML", package = "msdata")
mm8_sps <- Spectra(mm8_file)
mm8_db <- dbConnect(SQLite(), tempfile())
createMsBackendSqlDatabase(mm8_db, mm8_file, blob = FALSE)
mm8_be <- backendInitialize(MsBackendSql(), mm8_db)

mm8_db_blob <- dbConnect(SQLite(), tempfile())
createMsBackendSqlDatabase(mm8_db_blob, mm8_file, blob = TRUE)
mm8_be_blob <- backendInitialize(MsBackendSql(), mm8_db_blob)

mm14_file <- system.file("microtofq", "MM14.mzML", package = "msdata")
mm_db <- dbConnect(SQLite(), tempfile())
createMsBackendSqlDatabase(mm_db, c(mm8_file, mm14_file), blob = FALSE)
mm_be <- backendInitialize(MsBackendSql(), mm_db)

tmt_file <- proteomics(full.names = TRUE)[4L]
tmt_mzr <- backendInitialize(MsBackendMzR(), tmt_file)
tmt_db <- dbConnect(SQLite(), tempfile())
createMsBackendSqlDatabase(tmt_db, tmt_file, blob = FALSE)
tmt_be <- backendInitialize(MsBackendSql(), tmt_db)

test_check("MsBackendSql")

dbDisconnect(mm8_db)
