test_that(".valid_dbcon works", {
    expect_equal(.valid_dbcon(NULL), NULL)
    expect_match(.valid_dbcon(3), "expected to be")
    tf <- tempfile()
    cn <- dbConnect(SQLite(), tf)
    expect_match(.valid_dbcon(cn), "Database lacks")
    dbDisconnect(cn)
    unlink(tf)
})

test_that(".insert_data et al work", {
    db_file <- tempfile()
    db <- dbConnect(SQLite(), db_file)
    .insert_data(db, mm8_file, storage = "long")
    res <- dbListTables(db)
    expect_equal(res, c("msms_spectrum", "msms_spectrum_peak"))
    spd <- dbGetQuery(db, "select * from msms_spectrum")
    spd$smoothed <- as.logical(spd$smoothed)
    spd$centroided <- as.logical(spd$centroided)
    spd$spectrum_id_ <- NULL
    spd$dataStorage <- NULL
    ref <- as.data.frame(spectraData(mm8_sps))
    ref$dataStorage <- NULL
    expect_equal(ref, spd)

    pks <- dbGetQuery(db, "select * from msms_spectrum_peak")
    pks <- split(pks[, 1:2], pks$spectrum_id_)
    pks <- lapply(pks, function(z) {
        z <- as.matrix(z)
        rownames(z) <- NULL
        z
    })
    names(pks) <- NULL
    expect_equal(as.list(peaksData(mm8_sps)), pks)
    dbDisconnect(db)
    file.remove(db_file)

    db_file <- tempfile()
    db <- dbConnect(SQLite(), db_file)
    .insert_data(db, mm8_file, storage = "blob")
    res <- dbListTables(db)
    expect_equal(res, c("msms_spectrum", "msms_spectrum_peak_blob"))
    dbDisconnect(db)
    file.remove(db_file)

    db_file <- tempfile()
    db <- dbConnect(SQLite(), db_file)
    .insert_data(db, mm8_file, storage = "blob2")
    res <- dbListTables(db)
    expect_equal(res, c("msms_spectrum", "msms_spectrum_peak_blob2"))
    dbDisconnect(db)
    file.remove(db_file)

    ## Mocking a MySQL connection
    db_file <- tempfile()
    db <- dbConnect(SQLite(), db_file)
    db <- as(db, "DummySQL")
    with_mocked_bindings(
        ".is_maria_db" = function(x) TRUE,
        code = .insert_data(db, mm8_file)
    )
    dbDisconnect(db)
    file.remove(db_file)

    ## .insert_data mocking a duckdb
    tf <- tempfile()
    con_test <- dbConnect(SQLite(), tf)
    with_mocked_bindings(
        ".db_requires_peak_id" = function(x) TRUE,
        code = expect_true(
            .insert_data(con = con_test, c(mm8_file, mm14_file),
                         storage = "long"))
    )
    res <- dbGetQuery(con_test, "select * from msms_spectrum_peak")
    expect_equal(colnames(res), c("mz", "intensity", "spectrum_id_",
                                  "peak_id_"))
    expect_equal(res$peak_id_, seq_len(nrow(res)))
    dbDisconnect(con_test)
    unlink(tf)
})

test_that(".set_backend_insert_data works", {
    s <- Spectra(c(mm8_file, mm14_file))
    expect_error(.set_backend_insert_data(s, f = c(1, 2, 3)), "match length")

    ## blob2
    con_ref <- dbConnect(SQLite(), tempfile())
    createMsBackendSqlDatabase(con_ref, c(mm8_file, mm14_file), blob = TRUE,
                               peaksStorageMode = "blob2")
    be_ref <- backendInitialize(MsBackendSql(), dbcon = con_ref)

    con_test <- dbConnect(SQLite(), tempfile())
    .set_backend_insert_data(s, con = con_test, blob = TRUE,
                             peaksStorageMode = "blob2")
    be_test <- backendInitialize(MsBackendSql(), dbcon = con_test)
    expect_equal(length(be_ref), length(be_test))
    expect_equal(spectraData(be_ref, c("rtime", "dataOrigin")),
                 spectraData(be_test, c("rtime", "dataOrigin")))
    expect_equal(peaksData(be_ref), peaksData(be_test))
    expect_equal(be_ref$spectrum_id_, be_test$spectrum_id_)
    expect_equal(dbGetQuery(con_ref, "select * from msms_spectrum_peak_blob2"),
                 dbGetQuery(con_test, "select * from msms_spectrum_peak_blob2"))
    dbDisconnect(con_test)
    dbDisconnect(con_ref)

    ## blob; No chunk-wise processing
    con_ref <- dbConnect(SQLite(), tempfile())
    createMsBackendSqlDatabase(con_ref, c(mm8_file, mm14_file), blob = TRUE,
                               peaksStorageMode = "blob")
    be_ref <- backendInitialize(MsBackendSql(), dbcon = con_ref)
    con_test <- dbConnect(SQLite(), tempfile())
    .set_backend_insert_data(s, f = factor(), con = con_test,
                             peaksStorageMode = "blob")
    be_test <- backendInitialize(MsBackendSql(), dbcon = con_test)
    expect_equal(length(be_ref), length(be_test))
    expect_equal(spectraData(be_ref, c("rtime", "dataOrigin")),
                 spectraData(be_test, c("rtime", "dataOrigin")))
    expect_equal(peaksData(be_ref), peaksData(be_test))
    expect_equal(be_ref$spectrum_id_, be_test$spectrum_id_)
    expect_equal(dbGetQuery(con_ref, "select * from msms_spectrum_peak_blob"),
                 dbGetQuery(con_test, "select * from msms_spectrum_peak_blob"))
    dbDisconnect(con_test)
    dbDisconnect(con_ref)

    ## long; Arbitrary chunks.
    con_ref <- dbConnect(SQLite(), tempfile())
    createMsBackendSqlDatabase(con_ref, c(mm8_file, mm14_file), blob = FALSE)
    be_ref <- backendInitialize(MsBackendSql(), dbcon = con_ref)
    con_test <- dbConnect(SQLite(), tempfile())
    f <- sort(rep(1:10, length.out = length(s)))
    .set_backend_insert_data(s, f = f, con = con_test, blob = FALSE)
    be_test <- backendInitialize(MsBackendSql(), dbcon = con_test)
    expect_equal(length(be_ref), length(be_test))
    expect_equal(spectraData(be_ref, c("rtime", "dataOrigin")),
                 spectraData(be_test, c("rtime", "dataOrigin")))
    expect_equal(peaksData(be_ref), peaksData(be_test))
    expect_equal(be_ref$spectrum_id_, be_test$spectrum_id_)
    expect_equal(dbGetQuery(con_ref, "select * from msms_spectrum_peak"),
                 dbGetQuery(con_test, "select * from msms_spectrum_peak"))

    ## mock a MySQL connection
    dbDisconnect(con_test)
    con_test <- dbConnect(SQLite(), tempfile())
    con_test <- as(con_test, "DummySQL")
    with_mocked_bindings(
        ".insert_backend_blob2" = function(...) TRUE,
        ".is_maria_db" = function(x) TRUE,
        code = expect_true(.set_backend_insert_data(s, f = factor(),
                                                    con = con_test))
    )

    dbDisconnect(con_ref)
    dbDisconnect(con_test)

    ## mocking duckdb and adding peak_id_ works.
    tf <- tempfile()
    con_test <- dbConnect(SQLite(), tf)
    with_mocked_bindings(
        ".db_requires_peak_id" = function(x) TRUE,
        code = expect_true(
            .set_backend_insert_data(s, con = con_test, blob = FALSE,
                                     peaksStorageMode = "long"))
    )
    res <- dbGetQuery(con_test, "select * from msms_spectrum_peak")
    expect_equal(colnames(res), c("mz", "intensity", "spectrum_id_",
                                  "peak_id_"))
    expect_equal(res$peak_id_, seq_len(nrow(res)))
    dbDisconnect(con_test)
    unlink(tf)
})

test_that("createMsBackendSqlDatabase works", {
    cn <- dbConnect(SQLite(), tempfile())

    expect_false(createMsBackendSqlDatabase(cn))

    expect_error(createMsBackendSqlDatabase("b", "b"), "valid connection")
    expect_true(createMsBackendSqlDatabase(cn, mm8_file))

    dbDisconnect(cn)

    cn <- dbConnect(SQLite(), tempfile())
    expect_error(createMsBackendSqlDatabase(cn, "not existing"), "not found")

    dbDisconnect(cn)
})

test_that("MsBackendSql works", {
    res <- MsBackendSql()
    expect_s4_class(res, "MsBackendSql")
})

test_that(".fetch_peaks_data_long works", {
    expect_equal(.fetch_peaks_data_long(MsBackendSql()), list())
    expect_equal(
        .fetch_peaks_data_long(MsBackendSql(), columns = "intensity"), list())
    expect_equal(.fetch_peaks_data_long(MsBackendSql(), columns = "intensity",
                                        drop = TRUE), list())

    res <- .fetch_peaks_data_long(mm8_be_long)
    expect_true(is.list(res))
    ref <- peaksData(mm8_sps@backend)
    expect_equal(res, ref)
    res <- .fetch_peaks_data_long(mm8_be_long, columns = "intensity")
    expect_true(is.list(res))
    expect_true(is.matrix(res[[1]]))
    expect_equal(unlist(res), unlist(mm8_sps$intensity))
    res <- .fetch_peaks_data_long(mm8_be_long, columns = "intensity",
                                  drop = TRUE)
    expect_true(is.list(res))
    expect_equal(res, as.list(mm8_sps$intensity))

    ## subsets.
    a <- mm8_be_long[c(4, 1, 4, 5, 1, 3)]
    b <- mm8_sps@backend[c(4, 1, 4, 5, 1, 3)]
    expect_equal(.fetch_peaks_data_long(a), peaksData(b))
    expect_equal(.fetch_peaks_data_long(a, "mz", TRUE), as.list(b$mz))

    ## spectra without peaks data.
    a <- mm8_be_long[c(1:10)]
    a@spectraIds <- c(a@spectraIds, 200L)
    res <- .fetch_peaks_data_long(a)[[length(a) + 1L]]
    expect_true(is.matrix(res))
    expect_true(nrow(res) == 0L)
    expect_equal(colnames(res), c("mz", "intensity"))

    res <- .fetch_peaks_data_long(a, columns = "intensity")[[length(a) + 1L]]
    expect_true(is.matrix(res))
    expect_true(nrow(res) == 0L)
    expect_equal(colnames(res), c("intensity"))

    res <- .fetch_peaks_data_long(
        a, columns = "intensity", drop = TRUE)[[length(a) + 1L]]
    expect_true(is.numeric(res))
    expect_false(is.matrix(res))
    expect_equal(res, numeric())
})

test_that(".fetch_peaks_data_blob works", {
    expect_equal(.fetch_peaks_data_blob(MsBackendSql()), list())
    expect_equal(
        .fetch_peaks_data_blob(MsBackendSql(), columns = "intensity"), list())
    expect_equal(.fetch_peaks_data_blob(MsBackendSql(), columns = "intensity",
                                        drop = TRUE), list())

    res <- .fetch_peaks_data_blob(mm8_be_blob)
    expect_true(is.list(res))
    expect_true(is.matrix(res[[1L]]))
    expect_equal(colnames(res[[1L]]), c("mz", "intensity"))
    ref <- peaksData(mm8_sps@backend)
    expect_equal(res, ref)
    res <- .fetch_peaks_data_blob(mm8_be_blob,
                                  columns = c("intensity", "mz"))
    expect_true(is.list(res))
    expect_true(is.matrix(res[[1L]]))
    expect_equal(colnames(res[[1L]]), c("intensity", "mz"))

    res <- .fetch_peaks_data_blob(mm8_be_blob, columns = "intensity")
    expect_true(is.list(res))
    expect_true(is.matrix(res[[1L]]))
    expect_equal(colnames(res[[1L]]), "intensity")
    expect_equal(unlist(res), unlist(mm8_sps$intensity))
    res <- .fetch_peaks_data_blob(mm8_be_blob, columns = "intensity",
                                  drop = TRUE)
    expect_true(is.list(res))
    expect_equal(res, as.list(mm8_sps$intensity))

    ## subsets.
    a <- mm8_be_blob[c(4, 1, 4, 5, 1, 3)]
    b <- mm8_sps@backend[c(4, 1, 4, 5, 1, 3)]
    expect_equal(.fetch_peaks_data_blob(a), peaksData(b))
    expect_equal(.fetch_peaks_data_blob(a, "mz", TRUE), as.list(b$mz))
})

test_that(".fetch_peaks_data_blob2 works", {
    expect_equal(.fetch_peaks_data_blob2(MsBackendSql()), list())
    expect_equal(
        .fetch_peaks_data_blob2(MsBackendSql(), columns = "intensity"), list())
    expect_equal(.fetch_peaks_data_blob2(MsBackendSql(), columns = "intensity",
                                        drop = TRUE), list())

    res <- .fetch_peaks_data_blob2(mm8_be_blob2)
    expect_true(is.list(res))
    expect_true(is.matrix(res[[1L]]))
    expect_equal(colnames(res[[1L]]), c("mz", "intensity"))
    ref <- peaksData(mm8_sps@backend)
    expect_equal(res, ref)
    res <- .fetch_peaks_data_blob2(mm8_be_blob2,
                                   columns = c("intensity", "mz"))
    expect_true(is.list(res))
    expect_true(is.matrix(res[[1L]]))
    expect_equal(colnames(res[[1L]]), c("intensity", "mz"))

    res <- .fetch_peaks_data_blob2(mm8_be_blob2, columns = "intensity")
    expect_true(is.list(res))
    expect_true(is.matrix(res[[1L]]))
    expect_equal(colnames(res[[1L]]), "intensity")
    expect_equal(unlist(res), unlist(mm8_sps$intensity))
    res <- .fetch_peaks_data_blob2(mm8_be_blob2, columns = "intensity",
                                   drop = TRUE)
    expect_true(is.list(res))
    expect_equal(res, as.list(mm8_sps$intensity))

    ## subsets.
    a <- mm8_be_blob2[c(4, 1, 4, 5, 1, 3)]
    b <- mm8_sps@backend[c(4, 1, 4, 5, 1, 3)]
    expect_equal(.fetch_peaks_data_blob2(a), peaksData(b))
    expect_equal(.fetch_peaks_data_blob2(a, "mz", TRUE), as.list(b$mz))
})

test_that(".fetch_spectra_data_sql works", {
    res <- .fetch_spectra_data_sql(mm8_be_long, columns = c("rtime", "msLevel"))
    expect_true(is.data.frame(res))
    expect_identical(colnames(res), c("rtime", "msLevel"))
    expect_identical(length(mm8_be_long), nrow(res))
})

test_that(".disable_mysql_keys works", {
    ## Mocking the call since we don't have a MySQL database connection for
    ## testing
    expect_true(.disable_mysql_keys(new("DummySQL")))
})

test_that(".spectra_data_sql works", {
    expect_error(.spectra_data_sql(mm8_be_long, c("rtime", "other_col")),
                 "other_col not available.")
    res <- .spectra_data_sql(mm8_be_long, c("rtime", "msLevel", "mz"))
    expect_s4_class(res, "DataFrame")
    expect_identical(colnames(res), c("rtime", "msLevel", "mz"))
    expect_identical(length(mm8_be_long), nrow(res))
    expect_s4_class(res$mz, "NumericList")

    expect_error(
        .spectra_data_sql(mm8_be_long, columns = c("rtime", "other_column")),
        " not available.")

    tmp <- mm8_be_long[c(3, 2, 2, 3, 1, 10, 1)]
    res <- .spectra_data_sql(tmp, c("rtime", "msLevel", "mz"))
    expect_identical(colnames(res), c("rtime", "msLevel", "mz"))
    expect_identical(length(tmp), nrow(res))
    expect_s4_class(res$mz, "NumericList")

    tmp_sps <- mm8_sps[c(3, 2, 2, 3, 1, 10, 1)]
    expect_equal(tmp_sps$mz, res$mz)
    expect_equal(tmp_sps$rtime, res$rtime)
    expect_equal(tmp_sps$intensity, tmp$intensity)

    ## blob
    tmp <- mm8_be_blob[c(3, 2, 2, 3, 1, 10, 1)]
    res <- .spectra_data_sql(tmp, c("rtime", "msLevel", "mz"))
    expect_identical(colnames(res), c("rtime", "msLevel", "mz"))
    expect_identical(length(tmp), nrow(res))
    expect_s4_class(res$mz, "NumericList")

    expect_error(
        .spectra_data_sql(mm8_be_blob, columns = c("rtime", "other_column")),
        " not available.")

    expect_equal(tmp_sps$mz, res$mz)
    expect_equal(tmp_sps$rtime, res$rtime)
    expect_equal(tmp_sps$intensity, tmp$intensity)
})

test_that(".db_data_type works", {
    x <- data.frame(a = 1:4, b = TRUE, c = "TRUE")
    res <- .db_data_type(new("SQLiteConnection"), x)
    expect_equal(res, c(a = "INT", b = "SMALLINT", c = "TEXT"))
    setClass("MySQLConnection", contains = "SQLiteConnection")
    res <- .db_data_type(new("MySQLConnection"), x)
    expect_equal(res, c(a = "INTEGER", b = "INTEGER", c = "TEXT"))
})

test_that(".available_peaks_variables works", {
    res <- .available_peaks_variables(mm8_be_long)
    expect_equal(res, c("mz", "intensity"))

    res <- .available_peaks_variables(mm8_be_blob)
    expect_equal(res, c("mz", "intensity"))

    res <- .available_peaks_variables(mm8_be_blob2)
    expect_equal(res, c("mz", "intensity"))

    a <- mm8_be_blob2[integer()]
    res <- .available_peaks_variables(a)
    expect_equal(res, c("mz", "intensity"))

    res <- .available_peaks_variables(MsBackendSql())
    expect_equal(res, character())
})

test_that(".has_local_variable works", {
    res <- .has_local_variable(mm8_be_long, c("other_id"))
    expect_false(res)
    tmp <- mm8_be_long
    tmp$other_id <- "a"
    res <- .has_local_variable(tmp, c("other_id"))
    expect_true(res)
})

test_that(".is_maria_db works", {
    expect_false(.is_maria_db(10))
})

test_that(".is_duckdb works", {
    expect_false(.is_duckdb(10))
})

test_that(".blob_type works", {
    expect_equal(.blob_type("a"), "MEDIUMBLOB")
    with_mocked_bindings(
        ".is_duckdb" = function(x) TRUE,
        code = expect_equal(.blob_type("a"), "BLOB")
    )
})

test_that(".precursor_mz_query works", {
    res <- .precursor_mz_query(10, ppm = 0, tolerance = 0.1)
    expect_equal(res, "precursorMz >= 9.9 and precursorMz <= 10.1")

    res <- .precursor_mz_query(c(20, 10, 5), ppm = 0, tolerance = 0.1)
    expect_equal(res, paste0("precursorMz >= 19.9 and precursorMz <= 20.1 or ",
                             "precursorMz >= 9.9 and precursorMz <= 10.1 or ",
                             "precursorMz >= 4.9 and precursorMz <= 5.1"))

    res <- .precursor_mz_query(c(20, 10, 5), ppm = 100,
                               tolerance = c(0.1, 0.2, 0.1))
    expect_equal(res, paste0("precursorMz >= ", 20 - ppm(20, 100) - 0.1,
                             " and precursorMz <= ", 20 + ppm(20, 100) + 0.1,
                             " or precursorMz >= ", 10 - ppm(10, 100) - 0.2,
                             " and precursorMz <= ", 10 + ppm(10, 100) + 0.2,
                             " or precursorMz >= ", 5 - ppm(5, 100) - 0.1,
                             " and precursorMz <= ", 5 + ppm(5, 100) + 0.1))
})

test_that(".db_info_string works", {
    res <-.db_info_string(mm8_be_long)
    expect_true(is.character(res))
    expect_true(length(res) == 1L)
})

test_that(".combine works", {
    tmp <- split(mm8_be_long[1:10], 1:10)
    res <- .combine(tmp)
    expect_s4_class(res, "MsBackendSql")
    expect_equal(mm8_be_long[1:10], res)

    tmp <- .combine(list(mm8_be_long))
    expect_equal(length(tmp), length(mm8_be_long))
    expect_equal(rtime(tmp), rtime(mm8_be_long))

    a <- mm8_be_long[1:10]
    b <- setBackend(Spectra(a), MsBackendMemory())@backend
    expect_error(.combine(list(a, b)), "Can only merge backends of the same")

    expect_error(.combine(list(mm8_be_long, mm_be)), "connected to the same")
})

test_that(".initialize_tables works", {
    a <- new("DummySQL")
    .initialize_tables(a, cols = c(a = "TEXT"))
    with_mocked_bindings(
        ".is_maria_db" = function(x) TRUE,
        code = .initialize_tables(a, cols = c(a = "TEXT"))
    )

    tf <- tempfile()
    tmp <- dbConnect(SQLite(), tf)
    cols <- c(a = "INT", b = "DOUBLE", spectrumId = "TEXT")
    .initialize_tables(tmp, cols = cols)
    res <- dbGetQuery(tmp, "select * from msms_spectrum_peak")
    expect_equal(colnames(res), c("mz", "intensity", "spectrum_id_"))
    dbDisconnect(tmp)
    unlink(tf)
    ## peak_id_ column
    tf <- tempfile()
    tmp <- dbConnect(SQLite(), tf)
    with_mocked_bindings(
        ".db_requires_peak_id" = function(x) TRUE,
        code = .initialize_tables(tmp, cols = cols)
    )
    res <- dbGetQuery(tmp, "select * from msms_spectrum_peak")
    expect_equal(colnames(res), c("mz", "intensity", "spectrum_id_",
                                  "peak_id_"))
    dbDisconnect(tmp)
    unlink(tf)
})

test_that(".initialize_tables_sql works", {
    res <- .initialize_tables_sql(3, c("a", "b"))
    expect_true(length(res) == 2L)
    expect_equal(res[[2L]], paste0("CREATE TABLE msms_spectrum_peak (mz ",
                                   "DOUBLE, intensity REAL, spectrum_id_ ",
                                   "INTEGER);"))
    with_mocked_bindings(
        ".is_maria_db" = function(x) TRUE,
        code = {
            res <- .initialize_tables_sql(3, c("a", "b"))
        }
    )
    expect_match(res[[1L]], "ENGINE=ARIA;")
    expect_match(res[[2L]], "ENGINE=ARIA;")
    with_mocked_bindings(
        ".is_maria_db" = function(x) TRUE,
        code = {
            res <- .initialize_tables_sql(
                3, c("a", "b"), partitionBy = "spectrum")
        }
    )
    expect_match(res[[2L]], "ENGINE=ARIA PARTITION BY HASH (spectrum_id_",
                 fixed = TRUE)
    with_mocked_bindings(
        ".is_maria_db" = function(x) TRUE,
        code = {
            res <- .initialize_tables_sql(3, c("a", "b"), partitionBy = "chunk")
        }
    )
    expect_match(res[[2L]], "ENGINE=ARIA PARTITION BY HASH (partition_",
                 fixed = TRUE)
})

test_that(".load_data_file works", {
    d <- data.frame(a = 1:4, b = TRUE, c = FALSE, d = 5)
    with_mocked_bindings(
        "dbExecute" = function(...) TRUE,
        code = expect_true(.load_data_file(3, d))
    )
})

test_that(".insert_peaks works", {
    d <- data.frame(a = 1:4, b = TRUE, c = FALSE, d = 5)
    with_mocked_bindings(
        "dbExecute" = function(...) TRUE,
        ".is_maria_db" = function(x) TRUE,
        code = expect_true(.insert_peaks(3, d))
    )
})

test_that(".insert_spectra_variables works", {
    d <- data.frame(a = 1:4, b = TRUE, c = FALSE, d = 5)
    with_mocked_bindings(
        "dbExecute" = function(...) TRUE,
        ".is_maria_db" = function(x) TRUE,
        ".load_data_file" = function(con, data, name) {},
        code = expect_true(.insert_spectra_variables(mm_db, d))
    )
    s <- new("SQLiteConnection")
    with_mocked_bindings(
        "dbExecute" = function(...) TRUE,
        ".is_maria_db" = function(x) TRUE,
        ".load_data_file" = function(con, data, name) {},
        code = expect_true(.insert_spectra_variables(s, d))
    )
})

test_that(".insert_backend works", {
    with_mocked_bindings(
        ".insert_spectra_variables" = function(...) {},
        ".insert_peaks" = function(...) {},
        ".is_maria_db" = function(x) TRUE,
        code = expect_true(
            length(.insert_backend(3, mm8_sps, partitionBy = "chunk", 1L)) == 2L
        )
    )
    ## Insert in long form
    tf <- tempfile()
    tmp <- dbConnect(SQLite(), tf)
    sv <- spectraVariables(mm8_sps)
    spd <- as.data.frame(spectraData(mm8_sps[1:2], columns = sv))
    cols <- .db_data_type(tmp, spd)
    .initialize_tables(tmp, cols)
    res <- .insert_backend(tmp, mm8_sps, index = 77, peak_index = 13)
    expect_true(is.list(res))
    expect_equal(res$spectrum_id, length(mm8_sps) + 77)
    expect_equal(res$peak_id, 13)
    res <- dbGetQuery(tmp, "select * from msms_spectrum")
    expect_equal(res$spectrum_id_, seq(78, length.out = nrow(res)))
    res <- dbGetQuery(tmp, "select * from msms_spectrum_peak")
    expect_true(all(res$spectrum_id_ %in% seq(78, length.out = nrow(res))))
    dbDisconnect(tmp)
    unlink(tf)

    ## Insert with peak ID.
    tf <- tempfile()
    tmp <- dbConnect(SQLite(), tf)
    sv <- spectraVariables(mm8_sps)
    spd <- as.data.frame(spectraData(mm8_sps[1:2], columns = sv))
    cols <- .db_data_type(tmp, spd)
    with_mocked_bindings(
        ".db_requires_peak_id" = function(x) TRUE,
        code = {
            .initialize_tables(tmp, cols = cols)
            l <- .insert_backend(tmp, mm8_sps, index = 77, peak_index = 13)
        }
    )
    expect_true(is.list(l))
    expect_equal(l$spectrum_id, length(mm8_sps) + 77)
    res <- dbGetQuery(tmp, "select * from msms_spectrum")
    expect_equal(res$spectrum_id_, seq(78, length.out = nrow(res)))
    res <- dbGetQuery(tmp, "select * from msms_spectrum_peak")
    expect_true(all(res$spectrum_id_ %in% seq(78, length.out = nrow(res))))
    expect_equal(l$peak_id, nrow(res) + 13)
    expect_equal(res$peak_id_, seq(14, length.out  = nrow(res)))
    dbDisconnect(tmp)
    unlink(tf)
})


test_that(".initialize_tables_blob_sql works", {
    res <- .initialize_tables_blob_sql(3, c("a", "b"))
    expect_true(length(res) == 2L)
    expect_equal(res[[2L]], paste0("CREATE TABLE msms_spectrum_peak_blob (mz ",
                                   "MEDIUMBLOB, intensity MEDIUMBLOB, ",
                                   "spectrum_id_ INTEGER);"))
    with_mocked_bindings(
        ".is_maria_db" = function(x) TRUE,
        code = {
            res <- .initialize_tables_blob_sql(3, c("a", "b"))
        }
    )
    expect_match(res[[1L]], "ENGINE=ARIA;")
    expect_match(res[[2L]], "ENGINE=ARIA;")
    with_mocked_bindings(
        ".is_maria_db" = function(x) TRUE,
        code = {
            res <- .initialize_tables_blob_sql(3, c("a", "b"),
                                               partitionBy = "spectrum")
        }
    )
    expect_match(res[[2L]], "ENGINE=ARIA PARTITION BY HASH (spectrum_id_",
                 fixed = TRUE)
    with_mocked_bindings(
        ".is_maria_db" = function(x) TRUE,
        code = {
            res <- .initialize_tables_blob_sql(3, c("a", "b"),
                                               partitionBy = "chunk")
        }
    )
    expect_match(res[[2L]], "ENGINE=ARIA PARTITION BY HASH (partition_",
                 fixed = TRUE)
})

test_that(".create_from_spectra_data works", {
    ## wrong format or missing data.
    tmpf <- tempfile()
    tmpcon <- dbConnect(SQLite(), tmpf)
    dta <- spectraData(mm8_sps)
    expect_error(.create_from_spectra_data(tmpcon, dta), "required")

    ## blob
    dta <- spectraData(
        mm8_sps, columns = c(spectraVariables(mm8_sps), "mz", "intensity"))
    .create_from_spectra_data(tmpcon, dta, peaksStorageMode = "blob")
    res <- backendInitialize(MsBackendSql(), dbcon = tmpcon)
    expect_true(all(mm8_sps$dataStorage != res$dataStorage))
    expect_equal(rtime(mm8_be_long), rtime(res))
    expect_equal(mz(mm8_be_long), mz(res))
    expect_equal(intensity(mm8_be_long), intensity(res))
    tbls <- dbListTables(tmpcon)
    expect_equal(tbls, c("msms_spectrum", "msms_spectrum_peak_blob"))

    expect_error(.create_from_spectra_data(tmpcon, dta),
                 "contains already tables of a")
    dbDisconnect(tmpcon)

    tmpf <- tempfile()
    tmpcon <- dbConnect(SQLite(), tmpf)
    dta_2 <- dta[, !colnames(dta) %in% c("rtime", "msLevel")]
    .create_from_spectra_data(tmpcon, dta_2, blob = FALSE)
    res2 <- backendInitialize(MsBackendSql(), dbcon = tmpcon)
    expect_true(all(is.na(res2$msLevel)))
    expect_true(all(is.na(res2$rtime)))
    dbDisconnect(tmpcon)

    ## long format
    tmpf <- tempfile()
    tmpcon <- dbConnect(SQLite(), tmpf)
    .create_from_spectra_data(tmpcon, dta, blob = FALSE)
    tbls <- dbListTables(tmpcon)
    expect_equal(tbls, c("msms_spectrum", "msms_spectrum_peak"))
    res2 <- backendInitialize(MsBackendSql(), dbcon = tmpcon)
    expect_true(all(res2$dataStorage != mm8_be_long$dataStorage))
    expect_equal(rtime(res2), rtime(mm8_be_long))
    expect_equal(mz(res2), mz(mm8_be_long))
    expect_equal(intensity(res2), intensity(mm8_be_long))
    dbDisconnect(tmpcon)

    ## empty data frame
    tmpf <- tempfile()
    tmpcon <- dbConnect(SQLite(), tmpf)
    dta <- spectraData(
        mm8_sps[integer()],
        columns = c(spectraVariables(mm8_sps), "mz", "intensity"))
    .create_from_spectra_data(tmpcon, dta)
    res3 <- backendInitialize(MsBackendSql(), dbcon = tmpcon)
    expect_true(validObject(res3))
    expect_equal(spectraVariables(mm8_be_long), spectraVariables(res3))
    expect_equal(colnames(spectraData(mm8_be_long)),
                 colnames(spectraData(res3)))
    expect_true(length(res3) == 0L)
    dbDisconnect(tmpcon)

    ## spectra without m/z and intensity
    dta <- data.frame(msLevel = c(1L, 2L, 2L), rtime = c(12.2, 12.3, 13.1))
    tmpf <- tempfile()
    tmpcon <- dbConnect(SQLite(), tmpf)
    expect_error(.create_from_spectra_data(tmpcon, dta))
    dta$mz <- list(numeric(), numeric(), numeric())
    dta$intensity <- list(numeric(), numeric(), numeric())
    dta$my_col <- "a"
    .create_from_spectra_data(tmpcon, dta)
    res4 <- backendInitialize(MsBackendSql(), dbcon = tmpcon)
    expect_true(validObject(res4))
    expect_true(all(colnames(dta) %in% spectraVariables(res4)))
    tmp <- spectraData(res4)
    expect_equal(
        tmp$mz, IRanges::NumericList(list(numeric(), numeric(), numeric()),
                                     compress = FALSE))
    expect_equal(tmp$msLevel, dta$msLevel)
    expect_equal(tmp$rtime, dta$rtime)
    dbDisconnect(tmpcon)

    tmpf <- tempfile()
    tmpcon <- dbConnect(SQLite(), tmpf)
    tmpcon <- as(tmpcon, "DummySQL")
    dta <- spectraData(
        mm8_sps, columns = c(spectraVariables(mm8_sps), "mz", "intensity"))
    ## With mock to simulate MariaDB
    with_mocked_bindings(
        ".insert_spectra_variables" = function(...) TRUE,
        ".insert_peaks" = function(...) TRUE,
        ".is_maria_db" = function(x) TRUE,
        code = .create_from_spectra_data(tmpcon, dta, blob = FALSE)
    )
    dbDisconnect(tmpcon)
})

test_that(".drop_na_columns works", {
    tmp <- data.frame(a = 1:3, b = NA, d = 1:3)
    res <- .drop_na_columns(tmp)
    expect_equal(res, tmp[, c(1, 3)])

    tmp$a <- NA
    tmp$d <- NA
    res <- .drop_na_columns(tmp)
    expect_true(is.data.frame(tmp))
    expect_true(ncol(res) == 0)
    expect_true(nrow(res) == 3)

    tmp$a <- list(1:3, 2:3, 3)
    tmp$z <- "b"
    res <- .drop_na_columns(tmp)
    expect_equal(res, tmp[, c(1, 4)])

    ## want to keep specific columns with missing values.
    res <- .drop_na_columns(tmp, keep = c("d"))
    expect_equal(res, tmp[, c(1, 3, 4)])

    tmp$b <- "4"
    tmp$d <- 9
    res <- .drop_na_columns(tmp)
    expect_equal(res, tmp)
})

test_that("MsBackendSql works with duckdb", {
    ## only execute if duckdb is available.
    if (require("duckdb", quietly = TRUE)) {
        ## long format
        db_file <- tempfile()
        db <- dbConnect(duckdb(), db_file)
        .insert_data(db, mm8_file, storage = "long")
        res <- dbListTables(db)
        expect_equal(res, c("msms_spectrum", "msms_spectrum_peak"))
        tmp <- backendInitialize(MsBackendSql(), dbcon = db)
        expect_equal(rtime(tmp), rtime(mm8_sps))
        expect_equal(centroided(tmp), centroided(mm8_sps))
        expect_equal(msLevel(tmp), msLevel(mm8_sps))
        expect_equal(mz(tmp), mz(mm8_sps))
        expect_equal(intensity(tmp), intensity(mm8_sps))
        dbDisconnect(db)
        file.remove(db_file)

        ## blob2
        db_file <- tempfile()
        db <- dbConnect(duckdb(), db_file)
        .insert_data(db, mm8_file, storage = "blob2")
        res <- dbListTables(db)
        expect_equal(res, c("msms_spectrum", "msms_spectrum_peak_blob2"))
        tmp <- backendInitialize(MsBackendSql(), dbcon = db)
        expect_equal(rtime(tmp), rtime(mm8_sps))
        expect_equal(centroided(tmp), centroided(mm8_sps))
        expect_equal(msLevel(tmp), msLevel(mm8_sps))
        expect_equal(mz(tmp), mz(mm8_sps))
        expect_equal(intensity(tmp), intensity(mm8_sps))
        dbDisconnect(db)
        file.remove(db_file)
    }
})

test_that(".db_is_long_form works", {
    expect_equal(.db_is_long_form(MsBackendSql()), NA)
    expect_true(.db_is_long_form(mm8_be_long))
    expect_false(.db_is_long_form(mm8_be_blob))
    expect_false(.db_is_long_form(mm8_be_blob2))
})

test_that(".fetch_long_form_sql works", {
    tmp <- mm8_be_long
    ## only spectra variables
    res <- .fetch_long_form_sql(tmp, c("rtime", "msLevel", "scanIndex"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("rtime", "msLevel", "scanIndex"))
    expect_equal(res$rtime, tmp$rtime)
    expect_equal(res$scanIndex, tmp$scanIndex)
    ## arbitrary order
    idx <- c(4, 1, 9, 20)
    tmp <- mm8_be_long[idx]
    res <- .fetch_long_form_sql(tmp, c("rtime", "msLevel", "scanIndex"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("rtime", "msLevel", "scanIndex"))
    expect_equal(res$rtime, tmp$rtime)
    expect_equal(res$scanIndex, tmp$scanIndex)
    ## duplicated order
    idx <- c(4, 1, 4, 9, 20, 4)
    tmp <- mm8_be_long[idx]
    res <- .fetch_long_form_sql(tmp, c("rtime", "msLevel", "scanIndex"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("rtime", "msLevel", "scanIndex"))
    expect_equal(res$rtime, tmp$rtime)
    expect_equal(res$scanIndex, tmp$scanIndex)

    ## only peak variables
    tmp <- mm8_be_long
    res <- .fetch_long_form_sql(tmp, c("mz"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("mz"))
    expect_equal(res$mz, unlist(tmp$mz))

    ## arbitrary order
    idx <- c(3, 9, 1, 30, 5)
    tmp <- mm8_be_long[idx]
    res <- .fetch_long_form_sql(tmp, c("mz"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("mz"))
    expect_equal(res$mz, unlist(tmp$mz)) # Does not work

    ## duplicated order
    idx <- c(4, 1, 3, 1, 100, 1, 4)
    tmp <- mm8_be_long[idx]
    res <- .fetch_long_form_sql(tmp, c("mz"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("mz"))
    expect_equal(res$mz, unlist(tmp$mz))

    ## peak and spectra variables
    tmp <- mm8_be_long
    res <- .fetch_long_form_sql(tmp, c("scanIndex", "rtime", "intensity"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("scanIndex", "rtime", "intensity"))
    expect_equal(res$intensity, unlist(tmp$intensity))
    expect_equal(res$scanIndex, rep(tmp$scanIndex, lengths(tmp)))
    expect_equal(res$rtime, rep(tmp$rtime, lengths(tmp)))

    ## arbitrary order
    idx <- c(4, 9, 1, 3, 13)
    tmp <- mm8_be_long[idx]
    res <- .fetch_long_form_sql(tmp, c("scanIndex", "rtime", "intensity"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("scanIndex", "rtime", "intensity"))
    expect_equal(res$intensity, unlist(tmp$intensity))
    expect_equal(res$scanIndex, rep(tmp$scanIndex, lengths(tmp)))
    expect_equal(res$rtime, rep(tmp$rtime, lengths(tmp)))

    ## duplicated order
    idx <- c(4, 9, 1, 4, 13, 1, 9)
    tmp <- mm8_be_long[idx]
    res <- .fetch_long_form_sql(tmp, c("scanIndex", "rtime", "intensity"))
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("scanIndex", "rtime", "intensity"))
    expect_equal(res$intensity, unlist(tmp$intensity))
    expect_equal(res$scanIndex, rep(tmp$scanIndex, lengths(tmp)))
    expect_equal(res$rtime, rep(tmp$rtime, lengths(tmp)))

    ## Require peak_id_ LLLL mock .db_requires_peak_id...
})

test_that(".requires_peak_id works", {
    ## For not it returns only TRUE if a duckdb database is used.
    expect_false(.requires_peak_id(mm8_be_long))
    expect_false(.requires_peak_id(mm8_be_blob))
})

test_that(".has_peak_id works", {
    expect_false(.has_peak_id(mm8_be_blob))
    expect_false(.has_peak_id(mm8_be_long))
    expect_false(.has_peak_id(MsBackendSql()))
})

test_that(".db_requires_peak_id works", {
    expect_false(.db_requires_peak_id(mm8_be_long@dbcon))
})
