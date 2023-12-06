test_that(".insert_data et al work", {
    db_file <- tempfile()
    db <- dbConnect(SQLite(), db_file)
    .insert_data(db, mm8_file)
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
})

test_that(".set_backend_insert_data works", {
    s <- Spectra(c(mm8_file, mm14_file))
    expect_error(.set_backend_insert_data(object, f = c(1, 2, 3)))

    con_ref <- dbConnect(SQLite(), tempfile())
    createMsBackendSqlDatabase(con_ref, c(mm8_file, mm14_file))
    be_ref <- backendInitialize(MsBackendSql(), dbcon = con_ref)

    con_test <- dbConnect(SQLite(), tempfile())
    .set_backend_insert_data(s, con = con_test)
    be_test <- backendInitialize(MsBackendSql(), dbcon = con_test)
    expect_equal(length(be_ref), length(be_test))
    expect_equal(spectraData(be_ref, c("rtime", "dataOrigin")),
                 spectraData(be_test, c("rtime", "dataOrigin")))
    expect_equal(peaksData(be_ref), peaksData(be_test))
    expect_equal(be_ref$spectrum_id_, be_test$spectrum_id_)
    expect_equal(dbGetQuery(con_ref, "select * from msms_spectrum_peak_blob"),
                 dbGetQuery(con_test, "select * from msms_spectrum_peak_blob"))

    ## No chunk-wise processing
    dbDisconnect(con_test)
    con_test <- dbConnect(SQLite(), tempfile())
    .set_backend_insert_data(s, f = factor(), con = con_test)
    be_test <- backendInitialize(MsBackendSql(), dbcon = con_test)
    expect_equal(length(be_ref), length(be_test))
    expect_equal(spectraData(be_ref, c("rtime", "dataOrigin")),
                 spectraData(be_test, c("rtime", "dataOrigin")))
    expect_equal(peaksData(be_ref), peaksData(be_test))
    expect_equal(be_ref$spectrum_id_, be_test$spectrum_id_)
    expect_equal(dbGetQuery(con_ref, "select * from msms_spectrum_peak_blob"),
                 dbGetQuery(con_test, "select * from msms_spectrum_peak_blob"))

    ## Arbitrary chunks.
    dbDisconnect(con_test)
    con_test <- dbConnect(SQLite(), tempfile())
    f <- sort(rep(1:10, length.out = length(s)))
    .set_backend_insert_data(s, f = f, con = con_test)
    be_test <- backendInitialize(MsBackendSql(), dbcon = con_test)
    expect_equal(length(be_ref), length(be_test))
    expect_equal(spectraData(be_ref, c("rtime", "dataOrigin")),
                 spectraData(be_test, c("rtime", "dataOrigin")))
    expect_equal(peaksData(be_ref), peaksData(be_test))
    expect_equal(be_ref$spectrum_id_, be_test$spectrum_id_)
    expect_equal(dbGetQuery(con_ref, "select * from msms_spectrum_peak_blob"),
                 dbGetQuery(con_test, "select * from msms_spectrum_peak_blob"))

    dbDisconnect(con_ref)
    dbDisconnect(con_test)
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

test_that(".fetch_peaks_sql works", {
    res <- .fetch_peaks_sql(MsBackendSql(), columns = "intensity")
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 0)
    expect_identical(colnames(res), c("spectrum_id_", "intensity"))

    res <- .fetch_peaks_sql(mm8_be, columns = c("mz"))
    expect_true(is.data.frame(res))
    expect_identical(colnames(res), c("spectrum_id_", "mz"))
})

test_that(".fetch_peaks_sql_blob works", {
    res <- .fetch_peaks_sql_blob(MsBackendSql(), columns = "intensity")
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 0)
    expect_identical(colnames(res), c("spectrum_id_", "intensity"))

    res <- .fetch_peaks_sql_blob(mm8_be_blob, columns = "mz")
    expect_true(is.data.frame(res))
    expect_identical(colnames(res), c("spectrum_id_", "mz"))
    expect_true(is.list(res$mz))
})

test_that(".fetch_spectra_data_sql works", {
    res <- .fetch_spectra_data_sql(mm8_be, columns = c("rtime", "msLevel"))
    expect_true(is.data.frame(res))
    expect_identical(colnames(res), c("rtime", "msLevel"))
    expect_identical(length(mm8_be), nrow(res))
})

test_that(".spectra_data_sql works", {
    res <- .spectra_data_sql(mm8_be, c("rtime", "msLevel", "mz"))
    expect_s4_class(res, "DataFrame")
    expect_identical(colnames(res), c("rtime", "msLevel", "mz"))
    expect_identical(length(mm8_be), nrow(res))
    expect_s4_class(res$mz, "NumericList")

    tmp <- mm8_be[c(3, 2, 2, 3, 1, 10, 1)]
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

    expect_equal(tmp_sps$mz, res$mz)
    expect_equal(tmp_sps$rtime, res$rtime)
    expect_equal(tmp_sps$intensity, tmp$intensity)
})

test_that(".available_peaks_variables works", {
    res <- .available_peaks_variables(mm8_be)
    expect_equal(res, c("mz", "intensity"))
})

test_that(".has_local_variable works", {
    res <- .has_local_variable(mm8_be, c("other_id"))
    expect_false(res)
    tmp <- mm8_be
    tmp$other_id <- "a"
    res <- .has_local_variable(tmp, c("other_id"))
    expect_true(res)
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
    res <-.db_info_string(mm8_be)
    expect_true(is.character(res))
    expect_true(length(res) == 1L)
})

test_that(".combine works", {
    tmp <- split(mm8_be[1:10], 1:10)
    res <- .combine(tmp)
    expect_s4_class(res, "MsBackendSql")
    expect_equal(mm8_be[1:10], res)
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
    .create_from_spectra_data(tmpcon, dta)
    res <- backendInitialize(MsBackendSql(), dbcon = tmpcon)
    expect_true(all(mm8_sps$dataStorage != res$dataStorage))
    expect_equal(rtime(mm8_be), rtime(res))
    expect_equal(mz(mm8_be), mz(res))
    expect_equal(intensity(mm8_be), intensity(res))
    tbls <- dbListTables(tmpcon)
    expect_equal(tbls, c("msms_spectrum", "msms_spectrum_peak_blob"))

    ## long format
    tmpf <- tempfile()
    tmpcon <- dbConnect(SQLite(), tmpf)
    .create_from_spectra_data(tmpcon, dta, blob = FALSE)
    tbls <- dbListTables(tmpcon)
    expect_equal(tbls, c("msms_spectrum", "msms_spectrum_peak"))
    res2 <- backendInitialize(MsBackendSql(), dbcon = tmpcon)
    expect_true(all(res2$dataStorage != res$dataStorage))
    expect_equal(rtime(res2), rtime(res))
    expect_equal(mz(res2), mz(res))
    expect_equal(intensity(res2), intensity(res))

    ## empty data frame
    tmpf <- tempfile()
    tmpcon <- dbConnect(SQLite(), tmpf)
    dta <- spectraData(
        mm8_sps[integer()],
        columns = c(spectraVariables(mm8_sps), "mz", "intensity"))
    .create_from_spectra_data(tmpcon, dta)
    res3 <- backendInitialize(MsBackendSql(), dbcon = tmpcon)
    expect_true(validObject(res3))
    expect_equal(spectraVariables(res2), spectraVariables(res3))
    expect_equal(colnames(spectraData(res2)), colnames(spectraData(res3)))
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
