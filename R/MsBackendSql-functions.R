#' @rdname MsBackendSql
#'
#' @export MsBackendSql
MsBackendSql <- function() {
    new("MsBackendSql")
}

#' @importFrom DBI dbListTables
#'
#' @noRd
.valid_dbcon <- function(x) {
    if (length(x)) {
        if (!inherits(x, "DBIConnection"))
            return("'dbcon' is expected to be a connection to a database")
        tables <- dbListTables(x)
        if (!all(c("msms_spectrum", "msms_spectrum_peak") %in% tables |
                 c("msms_spectrum", "msms_spectrum_peak_blob") %in% tables))
            return("Database lacks some required tables.")
    }
    NULL
}

.dbcon <- function(x) {
    x@dbcon
}

#' Returns the spectra data, from the database and eventually filling with
#' *core* spectra variables, if they are not available in the database.
#'
#' The data can be either:
#' - in the database.
#' - in the local data (if new variables were added with $name <-).
#' - core spectra variables - if they are not in the database they have to be
#'   initialized with `NA` and the correct data type.
#'
#' @return a `data.frame` - always, even if only with a single column.
#'
#' @importFrom IRanges NumericList CharacterList
#'
#' @importFrom S4Vectors extractCOLS
#'
#' @importFrom S4Vectors make_zero_col_DFrame
#'
#' @importFrom methods as callNextMethod getMethod
#'
#' @author Johannes Rainer
#'
#' @noRd
.spectra_data_sql <- function(x, columns = spectraVariables(x)) {
    res <- getMethod("spectraData", "MsBackendCached")(x, columns = columns)
    if (is.null(res))
        res <- make_zero_col_DFrame(length(x))
    ## Define what needs to be still retrieved.
    db_cols <- intersect(columns, x@spectraVariables)
    db_cols <- db_cols[!db_cols %in% c("mz", "intensity", colnames(res))]
    mz_cols <- intersect(columns, c("mz", "intensity"))

    if (length(db_cols))
        res <- cbind(
            res, as(.fetch_spectra_data_sql(x, columns = db_cols), "DataFrame"))
    ## Get m/z and intensity values
    if (length(mz_cols)) {
        pks <-  x@peak_fun(x, columns = mz_cols)
        f <- factor(pks$spectrum_id_)
        if (any(mz_cols == "mz")) {
            if (is.numeric(pks$mz))
                mzs <- unname(split(pks$mz, f)[as.character(x@spectraIds)])
            else
                mzs <- pks$mz[match(x@spectraIds, pks$spectrum_id_)]
            res$mz <- NumericList(mzs, compress = FALSE)
        }
        if (any(mz_cols == "intensity")) {
            if (is.numeric(pks$intensity))
                ints <- unname(
                    split(pks$intensity, f)[as.character(x@spectraIds)])
            else
                ints <- pks$intensity[match(x@spectraIds, pks$spectrum_id_)]
            res$intensity <- NumericList(ints, compress = FALSE)
        }
    }
    if (!all(columns %in% colnames(res)))
        stop("Column(s) ", paste0(columns[!columns %in% names(res)],
                                  collapse = ", "), " not available.",
             call. = FALSE)
    if (any(columns == "centroided") && !is.logical(res$centroided))
        res$centroided <- as.logical(res$centroided)
    if (any(columns == "smoothed") && !is.logical(res$smoothed))
        res$smoothed <- as.logical(res$smoothed)
    extractCOLS(res, columns)
}

#' @importFrom DBI dbGetQuery
#'
#' @noRd
.fetch_peaks_sql <- function(x, columns = c("mz", "intensity")) {
    if (length(x@dbcon)) {
        dbGetQuery(
            x@dbcon,
            paste0("select spectrum_id_,", paste(columns, collapse = ","),
                   " from msms_spectrum_peak where spectrum_id_ in (",
                   paste0("'", unique(x@spectraIds), "'", collapse = ","),")"))
    } else {
        data.frame(spectrum_id_ = integer(), mz = numeric(),
                   intensity = numeric())[, c("spectrum_id_", columns)]
    }
}

.fetch_peaks_sql_blob <- function(x, columns = c("mz", "intensity")) {
    if (length(x@dbcon)) {
        res <- dbGetQuery(
            x@dbcon,
            paste0("select spectrum_id_,", paste(columns, collapse = ","),
                   " from msms_spectrum_peak_blob where spectrum_id_ in (",
                   paste0("'", unique(x@spectraIds), "'", collapse = ","),")"))
        if (any(colnames(res) == "mz"))
             res$mz <- lapply(res$mz, unserialize)
        if (any(colnames(res) == "intensity"))
            res$intensity <- lapply(res$intensity, unserialize)
        res
    } else {
        res <- data.frame(spectrum_id_ = integer())
        res$mz <- list()
        res$intensity <- list()
        res[, c("spectrum_id_", columns)]
    }
}

.fetch_spectra_data_sql <- function(x, columns = c("spectrum_id_")) {
    orig_columns <- columns
    sql_columns <- unique(c("spectrum_id_", columns))
    ## That turns out to be faster than dbBind if we use a field in the
    ## database that is unique (such as spectrum_id).
    res <- dbGetQuery(
        x@dbcon,
        paste0("select ", paste(sql_columns, collapse = ","), " from ",
               "msms_spectrum where spectrum_id_ in (",
               paste0("'", unique(x@spectraIds), "'", collapse = ", ") ,")"))
    idx <- match(x@spectraIds, res$spectrum_id_)
    res <- res[idx[!is.na(idx)], , drop = FALSE]
    rownames(res) <- NULL
    res[, orig_columns, drop = FALSE]
}

#' Get columns from the msms_spectrum_peak database table (dropping spectrum_id)
#'
#' @param x `MsBackendSql`
#'
#' @noRd
.available_peaks_variables <- function(x) {
    if (length(x@dbcon)) {
        tbl <- "msms_spectrum_peak"
        if (any(dbListTables(.dbcon(x)) == "msms_spectrum_peak_blob"))
            tbl <- "msms_spectrum_peak_blob"
        res <- dbGetQuery(
            .dbcon(x), paste0("select * from ", tbl, " limit 1"))
        colnames(res)[!colnames(res) %in% c("spectrum_id_", "peak_id")]
    } else character()
}

.is_maria_db <- function(x) {
    inherits(x, "MariaDBConnection")
}

##
## Insertion of data below.
##
.initialize_tables_sql <- function(con, cols, partitionBy = "none",
                                   partitionNumber = 10) {
    sql_a <- paste0("CREATE TABLE msms_spectrum (",
                    paste(names(cols), cols, collapse = ", "),
                    ", spectrum_id_ INTEGER, PRIMARY KEY (spectrum_id_))")
    sql_b <- paste0("CREATE TABLE msms_spectrum_peak (mz DOUBLE, intensity ",
                    "REAL, spectrum_id_ INTEGER")
    ## MySQL/MariaDB supports partitioning
    if (.is_maria_db(con)) {
        sql_a <- paste0(sql_a, " ENGINE=ARIA;")
        if (partitionBy == "none")
            sql_b <- paste0(sql_b, ", INDEX (spectrum_id_)) ENGINE=ARIA;")
        if (partitionBy == "spectrum")
            sql_b <- paste0(sql_b, ", INDEX (spectrum_id_)) ENGINE=ARIA ",
                            "PARTITION BY HASH (spectrum_id_) PARTITIONS ",
                            partitionNumber, ";")
        if (partitionBy == "chunk")
            sql_b <- paste0(sql_b, ", partition_ SMALLINT, ",
                            "INDEX (spectrum_id_)) ENGINE=ARIA ",
                            "PARTITION BY HASH (partition_) PARTITIONS ",
                            partitionNumber, ";")
    } else
        sql_b <- paste0(sql_b, ");")
    list(sql_a, sql_b)
}

.initialize_tables <- function(con, cols, partitionBy = "none",
                               partitionNumber = 10) {
    sql <- .initialize_tables_sql(con, cols, partitionBy, partitionNumber)
    res <- dbExecute(con, sql[[1L]])
    res <- dbExecute(con, sql[[2L]])
}

.initialize_tables_blob_sql <- function(con, cols, partitionBy = "none",
                                        partitionNumber = 10) {
    sql_a <- paste0("CREATE TABLE msms_spectrum (",
                    paste(names(cols), cols, collapse = ", "),
                    ", spectrum_id_ INTEGER, PRIMARY KEY (spectrum_id_))")
    sql_b <- paste0("CREATE TABLE msms_spectrum_peak_blob (mz MEDIUMBLOB, ",
                    "intensity MEDIUMBLOB, spectrum_id_ INTEGER")
    ## MySQL/MariaDB supports partitioning
    if (.is_maria_db(con)) {
        sql_a <- paste0(sql_a, " ENGINE=ARIA;")
        if (partitionBy == "none")
            sql_b <- paste0(sql_b, ", PRIMARY KEY (spectrum_id_)) ENGINE=ARIA;")
        if (partitionBy == "spectrum")
            sql_b <- paste0(sql_b, ", PRIMARY KEY (spectrum_id_)) ENGINE=ARIA ",
                            "PARTITION BY HASH (spectrum_id_) PARTITIONS ",
                            partitionNumber, ";")
        if (partitionBy == "chunk")
            sql_b <- paste0(sql_b, ", partition_ SMALLINT, ",
                            "PRIMARY KEY (spectrum_id_)) ENGINE=ARIA ",
                            "PARTITION BY HASH (partition_) PARTITIONS ",
                            partitionNumber, ";")
    } else
        sql_b <- paste0(sql_b, ");")
    list(sql_a, sql_b)
}

.initialize_tables_blob <- function(con, cols, partitionBy = "none",
                                    partitionNumber = 10) {
    sql <- .initialize_tables_blob_sql(con, cols, partitionBy, partitionNumber)
    res <- dbExecute(con, sql[[1L]])
    res <- dbExecute(con, sql[[2L]])
}

#' @importFrom DBI dbWriteTable
#'
#' @noRd
.insert_spectra_variables <- function(con, data) {
    info <- dbGetInfo(con)
    if (length(info$dbname))
        data$dataStorage <- info$dbname
    else data$dataStorage <- "<database>"
    if (.is_maria_db(con) || inherits(con, "MySQLConnection"))
        .load_data_file(con, data, "msms_spectrum")
    else
        dbWriteTable(con, name = "msms_spectrum", value = data, append = TRUE)
    invisible(TRUE)
}

.insert_peaks <- function(con, data) {
    if (.is_maria_db(con) || inherits(con, "MySQLConnection"))
        .load_data_file(con, data, "msms_spectrum_peak")
    else
        dbWriteTable(con, name = "msms_spectrum_peak",
                     value = data, append = TRUE)
    invisible(TRUE)
}

#' For MySQL databases: export data and use LOAD DATA FILE to import.
#'
#' @importFrom data.table fwrite
#'
#' @importFrom DBI dbExecute
#'
#' @noRd
.load_data_file <- function(con, data, name) {
    f <- tempfile()
    logicals <- which(vapply(data, is.logical, TRUE))
    if (length(logicals))
        for (i in logicals)
            data[, i] <- as.integer(data[, i])
    fwrite(data, file = f, row.names = FALSE, col.names = FALSE, sep = "\t",
           na = "\\N", eol = "\n", quote = FALSE, showProgress = FALSE)
    res <- dbExecute(
        con, paste0("LOAD DATA LOCAL INFILE '", f, "' INTO TABLE ", name,
                    " FIELDS TERMINATED BY 0x09;"))
    file.remove(f)
}

#' Inserts the data of a single backend to a database.
#'
#' @param con database
#'
#' @param x `MsBackend` or `Spectra`. Advantage of using a `Spectra` is that it
#'     performs by default parallel processing also on the `peaksData` call.
#'
#' @param index `integer(1)` defining the last used spectrum_id for peaks.
#'
#' @param partitionBy `character(1)` how and if the table should be partitioned.
#'
#' @param chunk `integer(1)` with the number of the current chunk.
#'
#' @return `integer(1)` last used spectrum_id
#'
#' @importFrom DBI dbDataType
#'
#' @noRd
.insert_backend <- function(con, x, index = 0L,
                            partitionBy = "none", chunk = 0L) {
    sv <- spectraVariables(x)
    sv <- sv[!sv %in% c("mz", "intensity")]
    spd <- as.data.frame(spectraData(x, columns = sv))
    spectrum_id <- seq(index + 1L, index + nrow(spd))
    spd$spectrum_id_ <- spectrum_id
    .insert_spectra_variables(con, spd)
    pks <- as(x, "list")
    lns <- lengths(pks) / 2
    pks <- as.data.frame(do.call(rbind, pks))
    pks$spectrum_id_ <- rep(spectrum_id, lns)
    if (partitionBy == "chunk" && .is_maria_db(con)) {
        ## Append an integer for the current processed chunk to be used for
        ## the partitioning
        pks$partition_ <- chunk
    }
    .insert_peaks(con, pks)
    spectrum_id[length(spectrum_id)]
}

#' store m/z and intensity values as BLOB.
#'
#' @noRd
.insert_backend_blob <- function(con, x, index = 0L, ...) {
    sv <- spectraVariables(x)
    sv <- sv[!sv %in% c("mz", "intensity")]
    spd <- as.data.frame(spectraData(x, columns = sv))
    spectrum_id <- seq(index + 1L, index + nrow(spd))
    spd$spectrum_id_ <- spectrum_id
    .insert_spectra_variables(con, spd)
    pks <- peaksData(x, columns = c("mz", "intensity"))
    template <- data.frame(matrix(ncol = 0, nrow = 1))
    pks_table <- do.call(rbind, lapply(pks, function(z) {
        template$mz <- list(serialize(z[, 1], NULL))
        template$intensity <- list(serialize(z[, 2], NULL))
        template
    }))
    pks_table$spectrum_id_ <- spectrum_id
    dbWriteTable(con, name = "msms_spectrum_peak_blob",
                 value = pks_table, append = TRUE)
    spectrum_id[length(spectrum_id)]
}

#' Insert data sequentially from files provided by x and adds it to a database
#'
#' @details
#'
#' MySQL databases support data partitioning hence splitting huge tables into
#' multiple *partitions* which can improve data insertion and index generation.
#' For very large datasets it is suggested to enable table partitioning by
#' selecting either `partitionBy = "spectrum"` or `partitionBy = "chunk"`. The
#' first assignes consecutive spectra to different partitions while the latter
#' puts spectra from files part of the same *chunk* into the same partition.
#' Both options have about the same performance but `partitionBy = "spectrum"`
#' requires less disk space.
#'
#' @param con database connection.
#'
#' @param partitionBy `character(1)` defining if and how the peak data table
#'     should be partitioned. `"none"`: no partitioning, `"spectrum"`: peaks
#'     are assigned to the partition based on the spectrum ID (number), i.e.
#'     spectra are evenly (consecutively) assigned across partitions. For
#'     `partitionNumber = 3`, the first spectrum is assigned to the first
#'     partition, the second to the second, the third to the third and the
#'     fourth spectrum again to the first partition. `"chunk"`:
#'     spectra processed as part of the same *chunk* are placed into the same
#'     partition. All spectra from the next processed chunk are assigned to the
#'     next partition. Note that this is only available for MySQL/MariaDB
#'     databases, i.e., if `con` is a `MariaDBConnection`.
#'     See details for more information.
#'
#' @param partitionNumber `integer(1)` defining the number of partitions the
#'     database table will be partitioned into (only supported for MySQL/MariaDB
#'     databases).
#'
#' @param chunksize `integer(1)` number of files that are imported/processed as
#'     one *chunk*. Smaller numbers will require less memory but might result in
#'     slower imports.
#'
#' @param x `character` with the raw data files from which the data should be
#'     imported.
#'
#' @param blob `logical(1)` whether m/z and intensity values should be stored
#'     as BLOB in the database.
#'
#' @importFrom progress progress_bar
#'
#' @importFrom BiocParallel bpparam
#'
#' @importFrom Spectra Spectra
#'
#' @noRd
.insert_data <- function(con, x, backend = MsBackendMzR(), chunksize = 10,
                         partitionBy = c("none", "spectrum", "chunk"),
                         partitionNumber = 10, blob = FALSE) {
    partitionBy <- match.arg(partitionBy)
    ## initialize backend and initialize database.
    be <- backendInitialize(backend, x[1L])
    sv <- spectraVariables(be)
    sv <- sv[!sv %in% c("mz", "intensity")]
    spd <- as.data.frame(spectraData(be, columns = sv))
    if (inherits(con, "MySQLConnection"))
        cols <- vapply(spd, function(z) dbDataType(con, z), character(1))
    else cols <- dbDataType(con, spd)
    if (blob) {
        .initialize_tables_blob(con, cols, partitionBy, partitionNumber)
        peak_table <- "msms_spectrum_peak_blob"
    } else {
        .initialize_tables(con, cols, partitionBy, partitionNumber)
        peak_table <- "msms_spectrum_peak"
    }
    ## Loop over x to insert data.
    index <- 0
    message("Importing data ... ")
    idxs <- seq_along(x)
    chunks <- split(idxs, ceiling(idxs / chunksize))
    pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                           "total (:percent) in ",
                                           ":elapsed"),
                           total = length(chunks), clear = FALSE, force = TRUE)
    pb$tick(0)
    if (.is_maria_db(con)) {
        res <- dbExecute(con, "SET FOREIGN_KEY_CHECKS = 0;")
        res <- dbExecute(con, "SET UNIQUE_CHECKS = 0;")
        res <- dbExecute(con, "ALTER TABLE msms_spectrum DISABLE KEYS;")
        res <- dbExecute(con,
                         paste0("ALTER TABLE ", peak_table, " DISABLE KEYS;"))
    }
    for (i in seq_along(chunks)) {
        s <- Spectra(source = backend, x[chunks[[i]]], BPPARAM = bpparam())
        if (blob)
            index <- .insert_backend_blob(con, s, index = index)
        else
            index <- .insert_backend(con, s, index = index, partitionBy, i)
        rm(s)
        gc()
        pb$tick(1)
    }
    .create_indices(con, peak_table)
}

#' Similar to the .insert_data but takes data from the provided `Spectra`
#' object and inserts that (chunk-wise) into the database.
#'
#' @noRd
.set_backend_insert_data <- function(object, f = processingChunkFactor(object),
                                     con, BPPARAM = SerialParam(),
                                     blob = TRUE, ...) {
    if (!length(f))
        f <- rep(1L, length(object))
    if (!is.factor(f))
        f <- force(factor(f, levels = unique(f)))
    if (length(f) != length(object))
        stop("length of 'f' has to match length of 'object'")
    sv <- spectraVariables(object)
    sv <- sv[!sv %in% c("mz", "intensity", "spectrum_id_")]
    spd <- as.data.frame(spectraData(object[1], columns = sv))
    if (inherits(con, "MySQLConnection"))
        cols <- vapply(spd, function(z) dbDataType(con, z), character(1))
    else cols <- dbDataType(con, spd)
    if (blob) {
        .initialize_tables_blob(con, cols, partitionBy = "none", 10)
        peak_table <- "msms_spectrum_peak_blob"
    } else {
        .initialize_tables(con, cols, partitionBy = "none", 10)
        peak_table <- "msms_spectrum_peak"
    }
    if (.is_maria_db(con)) {
        res <- dbExecute(con, "SET FOREIGN_KEY_CHECKS = 0;")
        res <- dbExecute(con, "SET UNIQUE_CHECKS = 0;")
        res <- dbExecute(con, "ALTER TABLE msms_spectrum DISABLE KEYS;")
        res <- dbExecute(con,
                         paste0("ALTER TABLE ", peak_table, " DISABLE KEYS;"))
    }
    index <- 0
    message("Importing data ... ")
    pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                           "total (:percent) in ",
                                           ":elapsed"),
                           total = length(levels(f)), clear = FALSE,
                           force = TRUE)
    pb$tick(0)
    for (l in levels(f)) {
        s <- Spectra(object@backend[f == l])
        if (blob)
            index <- .insert_backend_blob(con, s, index = index)
        else index <- .insert_backend(con, s, index = index)
        pb$tick(1)
    }
    .create_indices(con, peak_table)
}

.create_indices <- function(con, peak_table) {
    message("Creating indices ", appendLF = FALSE)
    if (.is_maria_db(con)) {
        res <- dbExecute(con, "SET FOREIGN_KEY_CHECKS = 1;")
        message(".", appendLF = FALSE)
        res <- dbExecute(con, "SET UNIQUE_CHECKS = 1;")
        message(".", appendLF = FALSE)
        res <- dbExecute(con, "ALTER TABLE msms_spectrum ENABLE KEYS;")
        message(".", appendLF = FALSE)
        res <- dbExecute(con, paste0("ALTER TABLE ",peak_table," ENABLE KEYS;"))
        message(".", appendLF = FALSE)
    } else {
        res <- dbExecute(con, paste0("CREATE INDEX peak_spectrum_id on ",
                                     peak_table, " (spectrum_id_)"))
        message(".", appendLF = FALSE)
        res <- dbExecute(con, paste0("CREATE INDEX spectrum_spectrum_id on ",
                                     "msms_spectrum (spectrum_id_)"))
        message(".", appendLF = FALSE)
    }
    ## create remaining indices
    res <- dbExecute(con, paste0("CREATE INDEX spectrum_rtime on ",
                                  "msms_spectrum (rtime)"))
    message(".", appendLF = FALSE)
    res <- dbExecute(con, paste0("CREATE INDEX spectrum_precursor_mz on ",
                                  "msms_spectrum (precursorMz)"))
    message(".", appendLF = FALSE)
    res <- dbExecute(con, paste0("CREATE INDEX spectrum_ms_level on ",
                                  "msms_spectrum (msLevel)"))
    message(" Done")
}

#' @rdname MsBackendSql
#'
#' @importFrom Spectra MsBackendMzR
#'
#' @export
createMsBackendSqlDatabase <- function(dbcon, x = character(),
                                       backend = MsBackendMzR(),
                                       chunksize = 10L, blob = TRUE,
                                       partitionBy = c("none", "spectrum",
                                                       "chunk"),
                                       partitionNumber = 10L) {
    partitionBy <- match.arg(partitionBy)
    if (!length(x)) return(FALSE)
    if (!inherits(dbcon, "DBIConnection"))
        stop("'dbcon' needs to be a valid connection to a database.")
    .insert_data(dbcon, x, backend, chunksize = chunksize,
                 partitionBy = partitionBy,
                 partitionNumber = partitionNumber[1L],
                 blob = blob)
    TRUE
}

.has_local_variable <- function(x, variable = character()) {
    all(variable %in% colnames(x@localData))
}

.subset_query <- function(object, qry) {
            ids <- dbGetQuery(.dbcon(object), qry)[, "spectrum_id_"]
            extractByIndex(object, which(object@spectraIds %in% ids))
}

.precursor_mz_query <- function(mz, ppm = 20, tolerance = 0) {
    lmz <- length(mz)
    if (length(ppm) != lmz)
        ppm <- rep(ppm[1L], lmz)
    if (length(tolerance) != lmz)
        tolerance <- rep(tolerance[1L], lmz)
    mzdiff <- ppm(mz, ppm) + tolerance
    mzr <- rep(mz, each = 2) + c(-1, 1) * rep(mzdiff, each = 2)
    qry <- paste0("precursorMz", c(" >= ", " <= "), mzr, c(" and ", " or "),
                  collapse = "")
    substring(qry, 1, nchar(qry) - 4)
}

#' Helper function to create the SQL query to fetch IDs
#'
#' @noRd
.id_query <- function(x) {
    qry <- paste0("select spectrum_id_ from msms_spectrum where ")
    if (length(x) < 2000000)
        qry <- paste0(qry, "spectrum_id_ in (",
                      paste0(x@spectraIds, collapse = ","), ") and ")
    qry
}

#' @importFrom MsCoreUtils vapply1c rbindFill
.combine <- function(objects) {
    if (length(objects) == 1L)
        return(objects[[1L]])
    if (!all(vapply1c(objects, class) == class(objects[[1]])))
        stop("Can only merge backends of the same type: ", class(objects[[1]]))
    ids <- vapply(objects, .db_info_string, character(1))
    if (length(unique(ids)) != 1L)
        stop("Can only merge backends connected to the same database")
    res <- objects[[1]]
    res@spectraIds <- unname(
        do.call(c, lapply(objects, function(z) z@spectraIds)))
    ## merge slots of MsBackendCached.
    res@localData <- do.call(
        rbindFill, lapply(objects, function(z) z@localData))
    if (!nrow(res@localData))
        res@localData <- data.frame(row.names = seq_along(res@spectraIds))
    res@spectraVariables <- unique(
        unlist(lapply(objects, function(z) z@spectraVariables),
               use.names = FALSE))
    res@nspectra <- length(res@spectraIds)
    res
}

#' @importFrom DBI dbGetInfo
.db_info_string <- function(x) {
    do.call(paste, dbGetInfo(.dbcon(x)))
}

#' Create a new `MsBackendSql` database and insert the (full) data provided
#' with a `spectraData` `DataFrame`. This is supposed to support changing
#' backend using `setBackend` to a `MsBackendSql`.
#'
#' @note
#'
#' At present we are not supporting additional peak variables. These
#' might eventually be provided with an additional parameter
#' `peaksVariables` and these would have to be extracted and
#' merged with the `peaks` data frame.
#'
#' @param data `DataFrame` with the full spectra data as returned by
#'     `spectraData`.
#'
#' @note
#'
#' There is some redundancy and code duplication between this function and
#' the `insert_backend` and `insert_backend_blob` functions. These functions
#' are not (yet) merged, because the latter are designed to be more memory
#' efficient and to enable import of large data set. This function assumes
#' that the full data is already loaded into memory and passed through the
#' `data` parameter to this function (`data` being a `data.frame` such as
#' returned by `spectraData`). The main reason why this function (or
#' `setBackend`) does not support parallel or chunk-wise insertion is the
#' primary keys of the spectra and the index creation.
#'
#' @author Johannes Rainer
#'
#' @importFrom Spectra coreSpectraVariables
#'
#' @noRd
.create_from_spectra_data <- function(dbcon, data, blob = TRUE, ...) {
    tbls <- dbListTables(dbcon)
    if (any(c("msms_spectrum", "msms_spectrum_peak",
              "msms_spectrum_peak_blob") %in% tbls))
        stop("'dbcon' contains already tables of a 'MsBackendSql' database. ",
             "If this error occurred during a 'setBackend' call, try ",
             "passing 'f = factor()' to that function.", call. = FALSE)
    if (!all(c("mz", "intensity") %in% colnames(data)))
        stop("'data' lacks required columns \"mz\" and \"intensity\"")
    if (any(colnames(data) == "spectrum_id_")) {
        warning("Replacing original column \"spectrum_id_\"")
        data$spectrum_id_ <- NULL
    }
    message("preparing data ... ", appendLF = FALSE)
    data <- .drop_na_columns(
        data, keep = setdiff(colnames(data), names(coreSpectraVariables())))
    if (!any(colnames(data) == "precursorMz"))
        data$precursorMz <- NA_real_
    if (!any(colnames(data) == "msLevel"))
        data$msLevel <- NA_integer_
    if (!any(colnames(data) == "rtime"))
        data$rtime <- NA_real_
    lns <- lengths(data$mz)
    mzs <- data$mz
    ints <- data$intensity
    data$mz <- NULL
    data$intensity <- NULL
    data <- as.data.frame(data)
    if (nrow(data))
        data$dataStorage <- "<database>"
    if (inherits(dbcon, "MySQLConnection"))
        cols <- vapply(data, function(z) dbDataType(dbcon, z), character(1))
    else cols <- dbDataType(dbcon, data)
    if (blob) {
        peak_table <- "msms_spectrum_peak_blob"
        .initialize_tables_blob(dbcon, cols)
    } else {
        peak_table <- "msms_spectrum_peak"
        .initialize_tables(dbcon, cols)
    }
    if (nrow(data)) {
        if (.is_maria_db(dbcon)) {
            res <- dbExecute(dbcon, "SET FOREIGN_KEY_CHECKS = 0;")
            res <- dbExecute(dbcon, "SET UNIQUE_CHECKS = 0;")
            res <- dbExecute(dbcon, "ALTER TABLE msms_spectrum DISABLE KEYS;")
            res <- dbExecute(
                dbcon, paste0("ALTER TABLE ", peak_table, " DISABLE KEYS;"))
        }
        sid <- seq_len(nrow(data))
        data$spectrum_id_ <- sid
        message("Done")
        message("Inserting data ... ", appendLF = FALSE)
        .insert_spectra_variables(dbcon, data)
        if (blob) {
            pks <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
            pks$mz <- lapply(as.list(mzs), serialize, connection = NULL)
            pks$intensity <- lapply(as.list(ints), serialize, connection = NULL)
            pks$spectrum_id_ <- sid
            dbWriteTable(dbcon, name = peak_table, value = pks, append = TRUE)
        } else {
            pks <- data.frame(mz = unlist(mzs), intensity = unlist(ints))
            pks$spectrum_id_ <- rep(sid, lns)
            .insert_peaks(dbcon, pks)
        }
        message("Done")
        .create_indices(dbcon, peak_table)
    } else message("Done \nInserting data ... Done")
}

#' @importFrom MsCoreUtils vapply1l
#'
#' @noRd
.drop_na_columns <- function(x, keep = character()) {
    if (!nrow(x))
        return(x)
    nas <- vapply1l(x, function(z) {
        allna <- all(is.na(z))
        if (length(allna) > 1)
            FALSE
        else allna
    })
    nas <- nas & !colnames(x) %in% keep
    if (any(nas))
        x[, !nas, drop = FALSE]
    else x
}
