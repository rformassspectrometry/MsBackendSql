test_that("backendInitialize works", {
    expect_error(backendInitialize(MsBackendSql()), "required")
    expect_error(backendInitialize(MsBackendSql(), dbcon = "file"), "connection")

    be <- backendInitialize(MsBackendSql(), dbcon = mm8_db)

    ## backendInitialize creating a new database.
    expect_warning(
        be2 <- backendInitialize(
            MsBackendSql(), dbcon = dbConnect(SQLite(), tempfile()),
            data = spectraData(be))
      , "Replacing")
    expect_equal(be$mz, be2$mz)
    expect_equal(be$rtime, be2$rtime)

    ## empty data
    expect_warning(
        be2 <- backendInitialize(
            MsBackendSql(), dbcon = dbConnect(SQLite(), tempfile()),
            data = spectraData(be[integer()]))
      , "Replacing")
    expect_true(length(be2) == 0L)
    expect_equal(spectraVariables(be2), spectraVariables(be))
})

test_that("dataStorage works", {
    res <- dataStorage(MsBackendSql())
    expect_identical(res, character())

    res <- dataStorage(mm8_be)
    expect_true(is.character(res))
    expect_identical(length(res), length(mm8_be))
})

test_that("extractByIndex,MsBackendSql works", {
    idx <- c(4L, 12L, 100L, 14L)
    res <- extractByIndex(mm8_be, idx)
    expect_identical(res@spectraIds, idx)
})

test_that("[,MsBackendSql works", {
    res <- mm8_be[]
    expect_equal(res, mm8_be)

    idx <- c(4L, 12L, 100L, 14L)
    res <- mm8_be[idx]
    expect_identical(res@spectraIds, idx)

    ## Duplicated elements.
    idx <- c(2L, 5L, 1L, 2L, 5L)
    res <- mm8_be[idx]
    expect_identical(res@spectraIds, idx)

    ## With additional data.
    be <- mm8_be
    be$new_var <- 1:length(be)
    res <- be[idx]
    expect_identical(be$new_var[idx], res$new_var)

    res <- peaksData(mm8_be[c(3, 1)])
    expect_true(length(res) == 2L)
    expect_equal(res[[1L]], peaksData(mm8_sps[3L])[[1L]])
    expect_equal(res[[2L]], peaksData(mm8_sps[1L])[[1L]])

    res <- peaksData(mm8_be_blob[c(3, 1)])
    expect_true(length(res) == 2L)
    expect_equal(res[[1L]], peaksData(mm8_sps[3L])[[1L]])
    expect_equal(res[[2L]], peaksData(mm8_sps[1L])[[1L]])

    res <- peaksData(mm8_be[1L])
    expect_true(length(res) == 1L)
    expect_true(is.matrix(res[[1L]]))
    expect_equal(res[[1L]], peaksData(mm8_sps[1L])[[1L]])

    res <- peaksData(mm8_be_blob[1])
    expect_true(length(res) == 1L)
    expect_true(is.matrix(res[[1L]]))
    expect_equal(res[[1L]], peaksData(mm8_sps[1L])[[1L]])
})

test_that("peaksData,MsBackendSql works", {
    idx <- c(4L, 12L, 100L, 14L)
    res <- mm8_be[idx]
    expect_identical(peaksData(res), peaksData(mm8_sps@backend[idx]))

    idx <- c(2L, 5L, 1L, 2L, 5L)
    res <- mm8_be[idx]
    expect_identical(peaksData(res), peaksData(mm8_sps@backend[idx]))

    res <- peaksData(mm8_be)
    res_2 <- peaksData(mm8_be, "mz")
    expect_true(colnames(res_2[[1]]) == "mz")
    expect_equal(res[[1]][, "mz"], res_2[[1]][, "mz"])

    res_2 <- peaksData(mm8_be, c("intensity", "mz"))
    expect_true(all(colnames(res_2[[1]]) == c("intensity", "mz")))
    expect_equal(res[[1]][, "mz"], res_2[[1]][, "mz"])
    expect_equal(res[[1]][, "intensity"], res_2[[1]][, "intensity"])

    ## blob
    idx <- c(4L, 12L, 100L, 14L)
    res <- mm8_be_blob[idx]
    expect_identical(peaksData(res), peaksData(mm8_sps@backend[idx]))

    idx <- c(2L, 5L, 1L, 2L, 5L)
    res <- mm8_be_blob[idx]
    expect_identical(peaksData(res), peaksData(mm8_sps@backend[idx]))

    res <- peaksData(mm8_be_blob)
    res_2 <- peaksData(mm8_be_blob, "mz")
    expect_true(colnames(res_2[[1]]) == "mz")
    expect_equal(res[[1]][, "mz"], res_2[[1]][, "mz"])

    res_2 <- peaksData(mm8_be_blob, c("intensity", "mz"))
    expect_true(all(colnames(res_2[[1]]) == c("intensity", "mz")))
    expect_equal(res[[1]][, "mz"], res_2[[1]][, "mz"])
    expect_equal(res[[1]][, "intensity"], res_2[[1]][, "intensity"])

    res_2 <- peaksData(mm8_be_blob[idx], c("intensity"))
    expect_equal(res_2, peaksData(mm8_sps@backend[idx], "intensity"))

    res <- peaksData(mm8_be[1L])
    expect_true(is.list(res))
    expect_true(length(res) == 1L)
    expect_true(is.matrix(res[[1L]]))
    expect_equal(colnames(res[[1L]]), c("mz", "intensity"))
})

test_that("peaksVariables,MsBackendSql works", {
    expect_equal(peaksVariables(mm8_be), c("mz", "intensity"))
})

test_that("intensity<-,MsBackendSql works", {
    expect_error(intensity(mm8_be) <- 1:5, "replace")
})

test_that("mz<-,MsBackendSql works", {
    expect_error(mz(mm8_be) <- 1:5, "replace")
})

test_that("spectraData,MsBackendSql works", {
    res <- spectraData(MsBackendSql())
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), spectraVariables(MsBackendSql()))

    res <- spectraData(MsBackendSql(), c("rtime", "mz"))
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), c("rtime", "mz"))

    res <- spectraData(mm8_be, c("msLevel", "rtime", "mz"))
    expect_equal(res, spectraData(mm8_sps, c("msLevel", "rtime", "mz")))

    ## Arbitrary ordering
    idx <- c(4L, 12L, 100L, 14L)
    be <- mm8_be[idx]
    expect_equal(spectraData(be, c("msLevel", "rtime", "mz")),
                 spectraData(mm8_sps@backend[idx], c("msLevel", "rtime", "mz")))

    idx <- c(2L, 5L, 1L, 2L, 5L)
    be <- mm8_be[idx]
    expect_equal(spectraData(be, c("msLevel", "rtime", "mz")),
                 spectraData(mm8_sps@backend[idx], c("msLevel", "rtime", "mz")))
})

test_that("$<-,MsBackendSql works", {
    be <- mm8_be
    expect_error(mm8_be$spectrum_id_ <- "a", "not be")
    be$new_var <- "A"
    expect_true(any(spectraVariables(be) == "new_var"))
    expect_true(all(be$new_var == "A"))
})

test_that("reset,MsBackendSql", {
    be <- mm8_be[c(5, 2, 10)]
    be$add_var <- "B"

    be_res <- reset(be)
    expect_identical(length(be_res), length(mm8_be))
})

test_that("spectraNames,spectraNames<-,MsBackendSql", {
    res <- spectraNames(mm8_be)
    expect_true(is.character(res))
    expect_identical(res, as.character(seq_along(mm8_be)))

    expect_error(spectraNames(mm8_be) <- rev(seq_along(mm8_be)),
                 "not supported")
})

test_that("filterMsLevel,MsBackendSql works", {
    res <- filterMsLevel(mm8_be)
    expect_equal(res, mm8_be)

    res <- filterMsLevel(mm8_be, msLevel = 1:2)
    expect_equal(res, mm8_be)

    res <- filterMsLevel(mm8_be, msLevel = 3)
    expect_true(length(res) == 0)

    tmp <- mm8_be
    tmp$msLevel <- rep(1:2, 99)
    res <- filterMsLevel(tmp, msLevel = 1L)
    expect_true(length(res) == (length(tmp) / 2))
})

test_that("filterRt,MsBackendSql works", {
    res <- filterRt(mm8_be)
    expect_equal(res, mm8_be)

    res <- filterRt(mm8_be, c(-Inf, Inf))
    expect_equal(res, mm8_be)

    res <- filterRt(mm8_be, c(-Inf, 50000))
    expect_true(length(res) == length(mm8_be))

    res <- filterRt(mm8_be, c(0, Inf))
    expect_true(length(res) == length(mm8_be))

    res <- filterRt(mm8_be, rt = c(1000, 2000))
    expect_true(length(res) == 0)

    res <- filterRt(mm8_be, rt = c(10, 20))
    expect_true(all(res$rtime > 10 & res$rtime < 20))

    res <- filterRt(mm8_be, rt = c(10, 20), msLevel. = 2)
    expect_equal(res, mm8_be)

    tmp <- mm8_be
    tmp$msLevel <- sample(1:3, length(tmp), replace = TRUE)
    res <- filterRt(tmp, rt = c(10, 20), msLevel. = 3)
    res_3 <- filterMsLevel(res, 3)
    expect_true(all(rtime(res_3) >= 10 & rtime(res_3) <= 20))
    expect_equal(filterMsLevel(res, c(1, 2)), filterMsLevel(tmp, c(1, 2)))

    ## TMT
    res <- filterRt(tmt_be, rt = c(200, 210), msLevel. = 2)
    res_2 <- filterMsLevel(res, 2)
    expect_true(all(rtime(res_2) >= 200 & rtime(res_2) <= 210))
    expect_false(all(rtime(res) >= 200 & rtime(res) <= 210))
    expect_equal(filterMsLevel(res, 1), filterMsLevel(tmt_be, 1))

    res2 <- filterRt(res, rt = c(205, 210), msLevel. = 1)
    res2_1 <- filterMsLevel(res2, 1)
    expect_true(all(rtime(res2_1) >= 205 & rtime(res2_1) <= 210))
    expect_equal(filterMsLevel(res2, 2), filterMsLevel(res, 2))
})

test_that("filterDataOrigin works", {
    res <- filterDataOrigin(mm_be, normalizePath(mm8_file))
    expect_true(all(res$dataOrigin == normalizePath(mm8_file)))

    res <- filterDataOrigin(mm_be, normalizePath(mm14_file))
    expect_true(all(res$dataOrigin == normalizePath(mm14_file)))

    res <- filterDataOrigin(mm_be, normalizePath(c(mm14_file, mm8_file)))
    expect_equal(unique(dataOrigin(res)), normalizePath(c(mm14_file, mm8_file)))
})

test_that("filterPrecursorMzRange works", {
    res <- filterPrecursorMzRange(tmt_be, c(660, 670))
    tmp <- tmt_be
    tmp$precursorMz <- precursorMz(tmt_be)
    res_2 <- filterPrecursorMzRange(tmp, c(660, 670))
    expect_equal(peaksData(res), peaksData(res_2))
    expect_true(all(precursorMz(res) >= 660 & precursorMz(res) <= 670))

    expect_equal(filterPrecursorMzRange(tmt_be), tmt_be)
})

test_that("filterPrecursorMzValues works", {
    tmt_be2 <- tmt_be
    tmt_be2$precursorMz <- tmt_be$precursorMz
    pmz <- c(620.1, 404.25, 417.7, 506.6)

    res <- filterPrecursorMzValues(tmt_be, pmz, tolerance = 0.1)
    res_2 <- filterPrecursorMzValues(tmt_be2, pmz, tolerance = 0.1)
    expect_equal(length(res), length(res_2))
    expect_equal(precursorMz(res), precursorMz(res_2))

    res_3 <- filterPrecursorMzValues(tmt_mzr, pmz, tolerance = 0.1)
    expect_equal(length(res), length(res_3))
    expect_equal(precursorMz(res), precursorMz(res_3))

    res <- filterPrecursorMzValues(tmt_be, sort(pmz), tolerance = 0.1)
    res_2 <- filterPrecursorMzValues(tmt_be2, sort(pmz), tolerance = 0.1)
    res_3 <- filterPrecursorMzValues(tmt_mzr, sort(pmz), tolerance = 0.1)
    expect_equal(precursorMz(res), precursorMz(res_2))
    expect_equal(precursorMz(res), precursorMz(res_3))

    pmz <- c(456.3, 503.7815)
    res <- filterPrecursorMzValues(tmt_be, pmz)
    res_2 <- filterPrecursorMzValues(tmt_be2, pmz)
    res_3 <- filterPrecursorMzValues(tmt_mzr, pmz)
    expect_equal(precursorMz(res), precursorMz(res_2))
    expect_equal(precursorMz(res), precursorMz(res_3))

    res <- filterPrecursorMzValues(tmt_be, pmz[c(2, 1)])
    res_2 <- filterPrecursorMzValues(tmt_be2, pmz[c(2, 1)])
    res_3 <- filterPrecursorMzValues(tmt_mzr, pmz[c(2, 1)])
    expect_equal(precursorMz(res), precursorMz(res_2))
    expect_equal(precursorMz(res), precursorMz(res_3))
})

test_that("uniqueMsLevels,MsBackendSql works", {
    expect_equal(uniqueMsLevels(tmt_be), unique(msLevel(tmt_be)))
    expect_equal(uniqueMsLevels(MsBackendSql()), integer())
})

test_that("backendMerge,MsBackendSql works", {
    empty <- mm8_be[integer()]
    res <- backendMerge(empty)
    expect_equal(res, empty)

    spl <- split(mm8_be[1:10], 1:10)
    spl[[5]] <- empty

    mm8_sub <- mm8_be[c(1, 2,3, 4, 6, 7, 8, 9, 10)]
    res <- backendMerge(spl)
    expect_s4_class(res, "MsBackendSql")
    expect_true(length(res) == 9L)
    expect_equal(rtime(res), rtime(mm8_sub))
    expect_equal(mz(res), mz(mm8_sub))

    spl[[2]]$other_var <- 2L
    res <- backendMerge(spl)
    expect_equal(res$other_var, c(NA, 2L, NA, NA, NA, NA, NA, NA, NA))
})

test_that("centroided,MsBackendSql works", {
    expect_true(is.logical(centroided(tmt_be)))
})

test_that("smoothed,MsBackendSql works", {
    expect_true(is.logical(smoothed(tmt_be)))
})

test_that("tic,MsBackendSql works", {
    res <- tic(tmt_be)
    expect_true(is.numeric(res))
    expect_true(all(!is.na(res)))

    res <- tic(mm_be)
    expect_true(is.numeric(res))
    expect_true(all(!is.na(res)))
    res_2 <- tic(mm_be, initial = FALSE)
    expect_true(sum(res != res_2) > 10)
})

test_that("supportsSetBackend,MsBackendSql works", {
    expect_true(supportsSetBackend(MsBackendSql()))
    expect_true(isReadOnly(MsBackendSql()))
})

test_that("setBackend works with MsBackendSql", {
    expect_error(setBackend(mm8_sps, MsBackendSql()), "required")
    tmpcon <- dbConnect(SQLite(), tempfile())
    expect_error(res <- setBackend(mm8_sps, MsBackendSql()), "dbcon")
    res <- setBackend(mm8_sps, MsBackendSql(), dbcon = tmpcon)
    expect_equal(dbListTables(tmpcon),
                 c("msms_spectrum", "msms_spectrum_peak_blob"))
    expect_equal(mz(res), mz(mm8_sps))
    expect_equal(rtime(res), rtime(mm8_sps))
    expect_s4_class(res@backend, "MsBackendSql")

    dbDisconnect(tmpcon)
    tmpcon <- dbConnect(SQLite(), tempfile())
    res <- setBackend(mm8_sps, MsBackendSql(),
                      dbcon = tmpcon, blob = FALSE)
    expect_equal(dbListTables(tmpcon),
                 c("msms_spectrum", "msms_spectrum_peak"))
    expect_equal(mz(res), mz(mm8_sps))
    expect_equal(rtime(res), rtime(mm8_sps))
    expect_s4_class(res@backend, "MsBackendSql")

    dbDisconnect(tmpcon)
    tmpcon <- dbConnect(SQLite(), tempfile())
    res <- setBackend(mm8_sps[integer()], MsBackendSql(),
                      dbcon = tmpcon, blob = FALSE)
    expect_equal(dbListTables(tmpcon),
                 c("msms_spectrum", "msms_spectrum_peak"))
    dbDisconnect(tmpcon)
})

test_that("backendBpparam,MsBackendSql works", {
    expect_s4_class(backendBpparam(MsBackendSql()), "SerialParam")
    expect_s4_class(backendBpparam(MsBackendSql(), MulticoreParam(2)),
                    "SerialParam")
})

test_that("setBackend,Spectra,MsBackendSql works", {
    ref <- Spectra(c(mm14_file, mm8_file))
    expect_error(setBackend(ref, MsBackendSql()), "'dbcon'")

    con_test <- dbConnect(SQLite(), tempfile())
    res <- setBackend(ref, MsBackendSql(), dbcon = con_test)
    expect_equal(spectraData(ref, c("rtime", "dataOrigin")),
                 spectraData(res, c("rtime", "dataOrigin")))
    expect_equal(peaksData(ref), peaksData(res))
    expect_true(length(processingLog(res)) > length(processingLog(ref)))

    dbDisconnect(con_test)
})
