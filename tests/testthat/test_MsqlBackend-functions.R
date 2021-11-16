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
