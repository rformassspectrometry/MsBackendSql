#' @rdname MsqlBackend
#'
#' @export MsqlBackend
MsqlBackend <- function() {
    new("MsqlBackend")
}

#' @importFrom DBI dbListTables
#'
#' @noRd
.valid_dbcon <- function(x) {
    if (length(x)) {
        if (!inherits(x, "DBIConnection"))
            return("'dbcon' is expected to be a connection to a database")
        tables <- dbListTables(x)
        if (!all(c("msms_spectrum", "msms_spectrum_peak") %in% tables))
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
        pks <- .fetch_peaks_sql(x, columns = mz_cols)
        f <- factor(pks$spectrum_id_)
        if (any(mz_cols == "mz")) {
            mzs <- unname(split(pks$mz, f)[as.character(x@spectraIds)])
            res$mz <- NumericList(mzs, compress = FALSE)
        }
        if (any(mz_cols == "intensity")) {
            ints <- unname(
                split(pks$intensity, f)[as.character(x@spectraIds)])
            res$intensity <- NumericList(ints, compress = FALSE)
        }
    }
    if (!all(columns %in% colnames(res)))
        stop("Column(s) ", paste0(columns[!columns %in% names(res)],
                                  collapse = ", "), " not available.",
             call. = FALSE)
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
#' @param x `MsqlBackend`
#'
#' @noRd
.available_peaks_variables <- function(x) {
    if (length(x@dbcon)) {
        res <- dbGetQuery(
            .dbcon(x), "select * from msms_spectrum_peak limit 1")
        colnames(res)[!colnames(res) %in% c("spectrum_id_", "peak_id")]
    } else character()
}

##
## Insertion of data below.
##

.initialize_tables <- function(con, cols, partitionBy = "none",
                               partitionNumber = 10) {
    sql_a <- paste0("CREATE TABLE msms_spectrum (",
                    paste(names(cols), cols, collapse = ", "),
                    ", spectrum_id_ INTEGER, PRIMARY KEY (spectrum_id_))")
    sql_b <- paste0("CREATE TABLE msms_spectrum_peak (mz DOUBLE, intensity ",
                    "REAL, spectrum_id_ INTEGER")
    ## MySQL/MariaDB supports partitioning
    if (inherits(con, "MariaDBConnection")) {
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
    suppressWarnings(res <- dbExecute(con, sql_a))
    suppressWarnings(res <- dbExecute(con, sql_b))
}

#' @importFrom DBI dbWriteTable
#'
#' @noRd
.insert_spectra_variables <- function(con, data) {
    if (inherits(con, "MariaDBConnection") || inherits(con, "MySQLConnection"))
        .load_data_file(con, data, "msms_spectrum")
    else
        dbWriteTable(con, name = "msms_spectrum", value = data, append = TRUE)
}

.insert_peaks <- function(con, data) {
    if (inherits(con, "MariaDBConnection") || inherits(con, "MySQLConnection"))
        .load_data_file(con, data, "msms_spectrum_peak")
    else
        dbWriteTable(con, name = "msms_spectrum_peak",
                     value = data, append = TRUE)
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
    logicals <- which(vapply(data, is.logical, logical(1)))
    if (length(logicals))
        for (i in logicals)
            data[, i] <- as.integer(data[, i])
    ## message("writing file")
    fwrite(data, file = f, row.names = FALSE, col.names = FALSE, sep = "\t",
           na = "\\N", eol = "\n", quote = FALSE, showProgress = FALSE)
    ## message("connection valid ", dbIsValid(conm))
    ## message("executing insert")
    res <- dbExecute(
        con, paste0("LOAD DATA LOCAL INFILE '", f, "' INTO TABLE ", name,
                    " FIELDS TERMINATED BY 0x09;"))
    ## message("removing file")
    res <- file.remove(f)
    if (!res)
        stop("failed to remove temporary file")
}

#' Inserts the data of a single backend to a database.
#'
#' @param con database
#'
#' @param x `MsBackend`
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
    pks <- peaksData(x)
    lns <- lengths(pks) / 2
    pks <- as.data.frame(do.call(rbind, pks))
    pks$spectrum_id_ <- rep(spectrum_id, lns)
    if (partitionBy == "chunk" && inherits(con, "MariaDBConnection")) {
        ## Append an integer for the current processed chunk to be used for
        ## the partitioning
        pks$partition_ <- chunk
    }
    .insert_peaks(con, pks)
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
#'     next partition. See details for more information.
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
#' @importFrom progress progress_bar
#'
#' @importFrom BiocParallel bpparam
#'
#' @noRd
.insert_data <- function(con, x, backend = MsBackendMzR(), chunksize = 10,
                         partitionBy = c("none", "spectrum", "chunk"),
                         partitionNumber = 10) {
    partitionBy <- match.arg(partitionBy)
    ## initialize backend and initialize database.
    be <- backendInitialize(backend, x[1L])
    sv <- spectraVariables(be)
    sv <- sv[!sv %in% c("mz", "intensity")]
    spd <- as.data.frame(spectraData(be, columns = sv))
    if (inherits(con, "MySQLConnection"))
        cols <- vapply(spd, function(z) dbDataType(con, z), character(1))
    else cols <- dbDataType(con, spd)
    .initialize_tables(con, cols, partitionBy, partitionNumber)
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
    if (inherits(con, "MariaDBConnection")) {
        res <- dbExecute(con, "SET FOREIGN_KEY_CHECKS = 0;")
        res <- dbExecute(con, "SET UNIQUE_CHECKS = 0;")
        res <- dbExecute(con, "ALTER TABLE msms_spectrum DISABLE KEYS;")
        res <- dbExecute(con, "ALTER TABLE msms_spectrum_peak DISABLE KEYS;")
    }
    for (i in seq_along(chunks)) {
        be <- backendInitialize(backend, x[chunks[[i]]], BPPARAM = bpparam())
        index <- .insert_backend(con, be, index = index, partitionBy, i)
        rm(be)
        gc()
        pb$tick(1)
    }
    message("Creating indices (this can take, depending on the database's size",
            ", very long) ", appendLF = FALSE)
    if (inherits(con, "MariaDBConnection")) {
        res <- dbExecute(con, "SET FOREIGN_KEY_CHECKS = 1;")
        message(".", appendLF = FALSE)
        res <- dbExecute(con, "SET UNIQUE_CHECKS = 1;")
        message(".", appendLF = FALSE)
        res <- dbExecute(con, "ALTER TABLE msms_spectrum ENABLE KEYS;")
        message(".", appendLF = FALSE)
        res <- dbExecute(con, "ALTER TABLE msms_spectrum_peak ENABLE KEYS;")
        message(".", appendLF = FALSE)
    } else {
        res <- dbExecute(con, paste0("CREATE INDEX peak_spectrum_id on ",
                                     "msms_spectrum_peak (spectrum_id_)"))
        message(".", appendLF = FALSE)
        res <- dbExecute(con, paste0("CREATE INDEX spectrum_spectrum_id on ",
                                     "msms_spectrum (spectrum_id_)"))
        message(".", appendLF = FALSE)
    }
    message(" Done")
}

#' @rdname MsqlBackend
#'
#' @importFrom Spectra MsBackendMzR
#'
#' @export
createMsqlBackendDatabase <- function(dbcon, x = character(),
                                      backend = MsBackendMzR(),
                                      chunksize = 10L) {
    if (!length(x)) return(FALSE)
    if (!inherits(dbcon, "DBIConnection"))
        stop("'dbcon' needs to be a valid connection to a database.")
    .insert_data(dbcon, x, backend, chunksize = chunksize,
                 partitionBy = "spectrum")
    TRUE
}

.has_local_variable <- function(x, variable = character()) {
    all(variable %in% colnames(x@localData))
}

.subset_query <- function(object, qry) {
            ids <- dbGetQuery(.dbcon(object), qry)[, "spectrum_id_"]
            object[object@spectraIds %in% ids]
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
