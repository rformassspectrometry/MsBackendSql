library(testthat)
library(MsBackendSql)
library(Spectra)
library(RSQLite)
library(MsDataHub)

setClass("DummySQL",
         contains = "SQLiteConnection")

setMethod("dbExecute", c("DummySQL", "character"),
          function(conn, statement, ...) {
              TRUE
          })

a_file <- MsDataHub::X20171016_POOL_POS_1_105.134.mzML()
a_sps <- Spectra(a_file)
a_db_long <- dbConnect(SQLite(), tempfile())
createMsBackendSqlDatabase(a_db_long, a_file, blob = FALSE,
                           peaksStorageMode = "long")
a_be_long <- backendInitialize(MsBackendSql(), a_db_long)

a_db_blob <- dbConnect(SQLite(), tempfile())
createMsBackendSqlDatabase(a_db_blob, a_file, blob = TRUE,
                           peaksStorageMode = "blob")
a_be_blob <- backendInitialize(MsBackendSql(), a_db_blob)

a_db_blob2 <- dbConnect(SQLite(), tempfile())
createMsBackendSqlDatabase(a_db_blob2, a_file, blob = TRUE,
                           peaksStorageMode = "blob2")
a_be_blob2 <- backendInitialize(MsBackendSql(), a_db_blob2)

################################################################################
##
##  OPTIONAL TESTS WITH duckdb
## library(duckdb)
## mm8_db_long <- dbConnect(duckdb(), tempfile())
## createMsBackendSqlDatabase(mm8_db_long, mm8_file, blob = FALSE,
##                            peaksStorageMode = "long")
## mm8_be_long <- backendInitialize(MsBackendSql(), mm8_db_long)

## mm8_db_blob <- dbConnect(duckdb(), tempfile())
## createMsBackendSqlDatabase(mm8_db_blob, mm8_file, blob = TRUE,
##                            peaksStorageMode = "blob")
## mm8_be_blob <- backendInitialize(MsBackendSql(), mm8_db_blob)

## mm8_db_blob2 <- dbConnect(duckdb(), tempfile())
## createMsBackendSqlDatabase(mm8_db_blob2, mm8_file, blob = TRUE,
##                            peaksStorageMode = "blob2")
## mm8_be_blob2 <- backendInitialize(MsBackendSql(), mm8_db_blob2)
##
################################################################################

b_file <- MsDataHub::X20171016_POOL_POS_3_105.134.mzML()
mm_db <- dbConnect(SQLite(), tempfile())
createMsBackendSqlDatabase(mm_db, c(a_file, b_file), blob = FALSE)
mm_be <- backendInitialize(MsBackendSql(), mm_db)

tmt_file <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()
tmt_mzr <- backendInitialize(MsBackendMzR(), tmt_file)
tmt_db <- dbConnect(SQLite(), tempfile())
createMsBackendSqlDatabase(tmt_db, tmt_file, blob = FALSE)
tmt_be <- backendInitialize(MsBackendSql(), tmt_db)

test_check("MsBackendSql")

test_suite <- system.file("test_backends", "test_MsBackend",
                          package = "Spectra")

be <- a_be_blob
test_dir(test_suite, stop_on_failure = TRUE)

be <- tmt_be[sample(seq_along(tmt_be), 300)]
test_dir(test_suite, stop_on_failure = TRUE)

dbDisconnect(a_db_long)
dbDisconnect(a_db_blob)
dbDisconnect(a_db_blob2)
dbDisconnect(mm_db)
dbDisconnect(tmt_db)
