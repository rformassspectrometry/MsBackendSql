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
