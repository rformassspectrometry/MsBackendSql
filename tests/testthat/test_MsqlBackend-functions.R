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
    expect_equal(as.data.frame(spectraData(mm8_sps)), spd)

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

test_that("createMsqlBackendDatabase works", {
    cn <- dbConnect(SQLite(), tempfile())

    expect_false(createMsqlBackendDatabase(cn))

    expect_error(createMsqlBackendDatabase("b", "b"), "valid connection")
    expect_true(createMsqlBackendDatabase(cn, mm8_file))

    dbDisconnect(cn)

    cn <- dbConnect(SQLite(), tempfile())
    expect_error(createMsqlBackendDatabase(cn, "not existing"), "not found")

    dbDisconnect(cn)
})

test_that("MsqlBackend works", {
    res <- MsqlBackend()
    expect_s4_class(res, "MsqlBackend")
})

test_that(".fetch_peaks_sql works", {
    res <- .fetch_peaks_sql(MsqlBackend(), columns = "intensity")
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 0)
    expect_identical(colnames(res), c("spectrum_id_", "intensity"))

    res <- .fetch_peaks_sql(mm8_be, columns = c("mz"))
    expect_true(is.data.frame(res))
    expect_identical(colnames(res), c("spectrum_id_", "mz"))
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
