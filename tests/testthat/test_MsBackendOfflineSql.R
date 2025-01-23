mm_be_off <- backendInitialize(MsBackendOfflineSql(), SQLite(),
                               dbname = dbGetInfo(mm_db)$dbname)

tmt_be_off <- backendInitialize(MsBackendOfflineSql(), SQLite(),
                                dbname = dbGetInfo(tmt_db)$dbname)

test_that("MsBackendOfflineSql works", {
    res <- MsBackendOfflineSql()
    expect_s4_class(res, "MsBackendOfflineSql")
    expect_true(validObject(res))
})

test_that(".db_connect works", {
    res <- .db_connect(MsBackendOfflineSql())
    expect_equal(res, NULL)
    res <- .db_connect(mm_be_off)
    expect_s4_class(res, "SQLiteConnection")
    dbDisconnect(res)
})

test_that("backendInitialize,MsBackendOfflineSql works", {
    expect_error(backendInitialize(MsBackendOfflineSql()), "must be specified")
    expect_error(backendInitialize(MsBackendOfflineSql(), SQLite()),
                 "At least the database name")

    dbn <- dbGetInfo(mm8_db_long)$dbname

    res <- backendInitialize(MsBackendOfflineSql(), SQLite(), dbname = dbn)
    expect_s4_class(res, "MsBackendOfflineSql")
    expect_false(dbIsValid(res@dbcon))
    expect_s4_class(res@driver, "SQLiteDriver")
    expect_equal(res@.tables, mm8_be_long@.tables)
    expect_equal(res@dbname, dbn)
    expect_true(validObject(res))

    expect_output(show(res), "MsBackendOfflineSql")

    ## with data. LLLLLL
    data <- spectraData(mm8_be_long)
    tf <- tempfile()
    res <- backendInitialize(MsBackendOfflineSql(),
                             drv = SQLite(), dbname = tf,
                             data = data)
    expect_s4_class(res, "MsBackendOfflineSql")
    expect_equal(rtime(res), data$rtime)
    unlink(tf)
})

test_that("dataStorage,MsBackendOfflineSql works", {
    expect_equal(basename(dataStorage(mm_be)), basename(dataStorage(mm_be_off)))
    expect_false(dbIsValid(mm_be_off@dbcon))
})

test_that("peaksData,MsBackendOfflineSql works", {
    res <- peaksData(mm_be_off)
    expect_false(dbIsValid(mm_be_off@dbcon))
    expect_equal(res, peaksData(mm_be))
})

test_that("[,MsBackendOfflineSql works", {
    idx <- c(4L, 12L, 100L, 14L)
    res <- mm_be_off[idx]
    expect_identical(res@spectraIds, idx)

    ## Duplicated elements.
    idx <- c(2L, 5L, 1L, 2L, 5L)
    res <- mm_be_off[idx]
    expect_identical(res@spectraIds, idx)

    ## With additional data.
    be <- mm_be_off
    expect_false(dbIsValid(be@dbcon))
    be$new_var <- 1:length(be)
    expect_false(dbIsValid(be@dbcon))

    res <- be[idx]
    expect_false(dbIsValid(res@dbcon))
    expect_identical(be$new_var[idx], res$new_var)

    res <- peaksData(mm_be_off[c(3, 1)])
    expect_true(length(res) == 2L)
    expect_equal(res[[1L]], peaksData(mm_be_off[3L])[[1L]])
    expect_equal(res[[2L]], peaksData(mm_be_off[1L])[[1L]])

    res <- peaksData(mm_be_off[1L])
    expect_true(length(res) == 1L)
    expect_true(is.matrix(res[[1L]]))
    expect_equal(res[[1L]], peaksData(mm_be[1L])[[1L]])
})

test_that("peaksData,MsBackendOfflineSql works", {
    idx <- c(4L, 12L, 100L, 14L)
    res <- mm_be_off[idx]
    expect_identical(peaksData(res), peaksData(mm_be[idx]))

    idx <- c(2L, 5L, 1L, 2L, 5L)
    res <- mm_be_off[idx]
    expect_identical(peaksData(res), peaksData(mm_be[idx]))

    res <- peaksData(mm_be_off)
    res_2 <- peaksData(mm_be_off, "mz")
    expect_true(colnames(res_2[[1]]) == "mz")
    expect_equal(res[[1]][, "mz"], res_2[[1]][, "mz"])

    res_2 <- peaksData(mm_be_off, c("intensity", "mz"))
    expect_true(all(colnames(res_2[[1]]) == c("intensity", "mz")))
    expect_equal(res[[1]][, "mz"], res_2[[1]][, "mz"])
    expect_equal(res[[1]][, "intensity"], res_2[[1]][, "intensity"])

    res <- peaksData(mm_be_off[1L])
    expect_true(is.list(res))
    expect_true(length(res) == 1L)
    expect_true(is.matrix(res[[1L]]))
    expect_equal(colnames(res[[1L]]), c("mz", "intensity"))
})

test_that("peaksVariables,MsBackendOfflineSql works", {
    expect_equal(peaksVariables(mm_be_off), c("mz", "intensity"))
    expect_false(dbIsValid(mm_be_off@dbcon))
})

test_that("intensity<-,MsBackendOfflineSql works", {
    expect_error(intensity(mm8_be_long) <- 1:5, "replace")
})

test_that("mz<-,MsBackendOfflineSql works", {
    expect_error(mz(mm8_be_long) <- 1:5, "replace")
})

test_that("spectraData,MsBackendOfflineSql works", {
    res <- spectraData(mm_be_off, c("msLevel", "rtime", "mz"))
    expect_equal(res, spectraData(mm_be, c("msLevel", "rtime", "mz")))

    ## Arbitrary ordering
    idx <- c(4L, 12L, 100L, 14L)
    be <- mm_be_off[idx]
    expect_equal(spectraData(be, c("msLevel", "rtime", "mz")),
                 spectraData(mm_be[idx], c("msLevel", "rtime", "mz")))

    idx <- c(2L, 5L, 1L, 2L, 5L)
    be <- mm_be_off[idx]
    expect_equal(spectraData(be, c("msLevel", "rtime", "mz")),
                 spectraData(mm_be[idx], c("msLevel", "rtime", "mz")))
})

test_that("$<-,MsBackendOfflineSql works", {
    be <- mm_be_off
    expect_error(be$spectrum_id_ <- "a", "not be")
    be$new_var <- "A"
    expect_false(dbIsValid(be@dbcon))
    expect_true(any(spectraVariables(be) == "new_var"))
    expect_true(all(be$new_var == "A"))
    expect_false(dbIsValid(be@dbcon))
})

test_that("reset,MsBackendOfflineSql", {
    be <- mm_be_off[c(5, 2, 10)]
    be$add_var <- "B"

    be_res <- reset(be)
    expect_identical(length(be_res), length(mm_be))
})

test_that("spectraNames,spectraNames<-,MsBackendOfflineSql", {
    res <- spectraNames(mm_be_off)
    expect_true(is.character(res))
    expect_identical(res, as.character(seq_along(mm_be)))

    expect_error(spectraNames(mm_be_off) <- rev(seq_along(mm_be)),
                 "not supported")
})

test_that("filterMsLevel,MsBackendOfflineSql works", {
    res <- filterMsLevel(mm_be_off)
    expect_equal(res, mm_be_off)
    expect_false(dbIsValid(res@dbcon))

    res <- filterMsLevel(mm_be_off, msLevel = 1:2)
    expect_equal(res, mm_be_off)

    res <- filterMsLevel(mm_be_off, msLevel = 3)
    expect_true(length(res) == 0)

    tmp <- mm_be_off
    tmp$msLevel <- rep(1:2, length(tmp)/2)
    res <- filterMsLevel(tmp, msLevel = 1L)
    expect_true(length(res) == (length(tmp) / 2))
    expect_false(dbIsValid(res@dbcon))
})

test_that("filterRt,MsBackendOfflineSql works", {
    res <- filterRt(mm_be_off)
    expect_equal(res, mm_be_off)

    res <- filterRt(mm_be_off, rt = c(1000, 2000))
    expect_true(length(res) == 0)
    expect_false(dbIsValid(res@dbcon))

    res <- filterRt(mm_be_off, rt = c(10, 20))
    expect_true(all(res$rtime > 10 & res$rtime < 20))
    expect_false(dbIsValid(res@dbcon))

    res <- filterRt(mm_be_off, rt = c(10, 20), msLevel. = 2)
    expect_equal(res, mm_be_off)
    expect_false(dbIsValid(res@dbcon))

    tmp <- mm_be_off
    tmp$msLevel <- sample(1:3, length(tmp), replace = TRUE)
    res <- filterRt(tmp, rt = c(10, 20), msLevel. = 3)
    res_3 <- filterMsLevel(res, 3)
    expect_true(all(rtime(res_3) >= 10 & rtime(res_3) <= 20))
    expect_equal(filterMsLevel(res, c(1, 2)), filterMsLevel(tmp, c(1, 2)))
    expect_false(dbIsValid(res@dbcon))
})

test_that("filterDataOrigin,MsBackendOfflineSql works", {
    res <- filterDataOrigin(mm_be_off, normalizePath(mm8_file))
    expect_true(all(res$dataOrigin == normalizePath(mm8_file)))

    res <- filterDataOrigin(mm_be_off, normalizePath(mm14_file))
    expect_true(all(res$dataOrigin == normalizePath(mm14_file)))

    res <- filterDataOrigin(mm_be_off, normalizePath(c(mm14_file, mm8_file)))
    expect_equal(unique(dataOrigin(res)), normalizePath(c(mm14_file, mm8_file)))
})

test_that("filterPrecursorMzRange,MsBackendOfflineSql works", {
    res <- filterPrecursorMzRange(mm_be_off, c(100, 200))
    expect_s4_class(res, "MsBackendOfflineSql")
    expect_true(length(res) == 0)

    res <- filterPrecursorMzRange(tmt_be_off, c(500, 600))
    expect_s4_class(res, "MsBackendOfflineSql")
    expect_true(length(res) > 0)
    expect_true(all(msLevel(res) == 2L))
    expect_true(all(precursorMz(res) > 500 & precursorMz(res) < 600))
})

test_that("filterPrecursorMzValues,MsBackendOfflineSql works", {
    res <- filterPrecursorMzValues(tmt_be_off, 517, tolerance = 1)
    expect_s4_class(res, "MsBackendOfflineSql")
    expect_true(length(res) > 0)
    expect_true(all(msLevel(res) == 2L))
    expect_true(all(precursorMz(res) > 515 & precursorMz(res) < 519))
})

test_that("uniqueMsLevels,MsBackendOfflineSql works", {
    expect_equal(uniqueMsLevels(mm_be_off), 1L)
    expect_false(dbIsValid(mm_be_off@dbcon))
})

test_that("backendMerge,MsBackendOfflineSql works", {
    empty <- mm_be_off[integer()]
    res <- backendMerge(empty)
    expect_equal(res, empty)

    spl <- split(mm_be_off[1:10], 1:10)
    spl[[5]] <- empty

    mm_be_off_sub <- mm_be_off[c(1, 2,3, 4, 6, 7, 8, 9, 10)]
    res <- backendMerge(spl)
    expect_false(dbIsValid(res@dbcon))
    expect_s4_class(res, "MsBackendOfflineSql")
    expect_true(length(res) == 9L)
    expect_equal(rtime(res), rtime(mm_be_off_sub))
    expect_equal(mz(res), mz(mm_be_off_sub))

    spl[[2]]$other_var <- 2L
    res <- backendMerge(spl)
    expect_equal(res$other_var, c(NA, 2L, NA, NA, NA, NA, NA, NA, NA))
})

test_that("centroided,MsBackendOfflineSql works", {
    expect_true(is.logical(centroided(mm_be_off)))
})

test_that("smoothed,MsBackendOfflineSql works", {
    expect_true(is.logical(smoothed(mm_be_off)))
})

test_that("tic,MsBackendOfflineSql works", {
    res <- tic(mm_be_off)
    expect_true(is.numeric(res))
    expect_true(all(!is.na(res)))

    expect_true(all(!is.na(res)))
    res_2 <- tic(mm_be_off, initial = FALSE)
    expect_true(sum(res != res_2) > 10)
})

test_that("supportsSetBackend,MsBackendOfflineSql works", {
    expect_true(supportsSetBackend(MsBackendOfflineSql()))
    expect_true(isReadOnly(MsBackendOfflineSql()))
})

test_that("setBackend works with MsBackendOfflineSql", {
    expect_error(setBackend(mm8_sps, MsBackendOfflineSql()), "'drv'")
    expect_error(setBackend(mm8_sps, MsBackendOfflineSql(), drv = SQLite()),
                 "'dbname'")
    dbn <- tempfile()
    res <- setBackend(mm8_sps, MsBackendOfflineSql(), drv = SQLite(),
                      dbname = dbn)
    expect_true(inherits(res@backend, "MsBackendOfflineSql"))
    expect_equal(rtime(res), rtime(mm8_sps))
    expect_equal(peaksData(res), peaksData(mm8_sps))
})

test_that("backendBpparam,MsBackendOfflineSql works", {
    mcp <- MulticoreParam(2)
    expect_equal(backendBpparam(MsBackendOfflineSql(), mcp), mcp)
})

test_that("setBackend,Spectra,MsBackendOfflineSql works", {
    ref <- Spectra(c(mm14_file, mm8_file))
    expect_error(setBackend(ref, MsBackendOfflineSql()), "'drv'")
    expect_error(setBackend(ref, MsBackendOfflineSql(), drv = SQLite()),
                 "'dbname'")

    dbname_test <- tempfile()
    res <- setBackend(ref, MsBackendOfflineSql(), drv = SQLite(),
                      dbname = dbname_test)
    expect_equal(spectraData(ref, c("rtime", "dataOrigin")),
                 spectraData(res, c("rtime", "dataOrigin")))
    expect_equal(peaksData(ref), peaksData(res))
    expect_true(length(processingLog(res)) > length(processingLog(ref)))

    ref <- Spectra()
    dbf <- tempfile()
    res <- setBackend(ref, MsBackendOfflineSql(), drv = SQLite(), dbname = dbf)
    expect_s4_class(res, "Spectra")
    expect_true(length(res) == 0)
    expect_equal(msLevel(ref), msLevel(res))
    unlink(dbf)
})

test_that("filterRt,MsBackendOfflineSql works properly", {
    ref <- Spectra(c(mm14_file, mm8_file))
    dbname_test <- tempfile()
    res <- setBackend(ref, MsBackendOfflineSql(), drv = SQLite(),
                      dbname = dbname_test)
    expect_output(show(ref), "MsBackendMzR")
    expect_output(show(res), "MsBackendOfflineSql")
})
