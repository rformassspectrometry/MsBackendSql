test_that("backendInitialize works", {
    expect_error(backendInitialize(MsqlBackend()), "required")
    expect_error(backendInitialize(MsqlBackend(), dbcon = "file"), "connection")

    be <- backendInitialize(MsqlBackend(), dbcon = mm8_db)
})

test_that("dataStorage works", {
    res <- dataStorage(MsqlBackend())
    expect_identical(res, character())

    res <- dataStorage(mm8_be)
    expect_true(is.character(res))
    expect_identical(length(res), length(mm8_be))
})

test_that("[,MsqlBackend works", {
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
})

test_that("peaksData,MsqlBackend works", {
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
})

test_that("peaksVariables,MsqlBackend works", {
    expect_equal(peaksVariables(mm8_be), c("mz", "intensity"))
})

test_that("intensity<-,MsqlBackend works", {
    expect_error(intensity(mm8_be) <- 1:5, "replace")
})

test_that("mz<-,MsqlBackend works", {
    expect_error(mz(mm8_be) <- 1:5, "replace")
})

test_that("spectraData,MsqlBackend works", {
    res <- spectraData(MsqlBackend())
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), spectraVariables(MsqlBackend()))

    res <- spectraData(MsqlBackend(), c("rtime", "mz"))
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

test_that("$<-,MsqlBackend works", {
    be <- mm8_be
    expect_error(mm8_be$spectrum_id_ <- "a", "not be")
    be$new_var <- "A"
    expect_true(any(spectraVariables(be) == "new_var"))
    expect_true(all(be$new_var == "A"))
})

test_that("reset,MsqlBackend", {
    be <- mm8_be[c(5, 2, 10)]
    be$add_var <- "B"

    be_res <- reset(be)
    expect_identical(length(be_res), length(mm8_be))
})

test_that("spectraNames,spectraNames<-,MsqlBackend", {
    res <- spectraNames(mm8_be)
    expect_true(is.character(res))
    expect_identical(res, as.character(seq_along(mm8_be)))

    expect_error(spectraNames(mm8_be) <- 1:20, "does not support")
})

test_that("filterMsLevel,MsqlBackend works", {
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

test_that("filterRt,MsqlBackend works", {
    res <- filterRt(mm8_be)
    expect_equal(res, mm8_be)

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
