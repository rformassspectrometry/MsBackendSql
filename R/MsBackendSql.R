#' @include MsBackendSql-functions.R

#' @title `Spectra` MS backend storing data in a SQL database
#'
#' @aliases MsBackendSql-class
#' @aliases backendMerge,MsBackendOfflineSql-method
#' @aliases dataStorage,MsBackendOfflineSql-method
#' @aliases dbconn,MsBackendOfflineSql-method
#' @aliases filterDataOrigin,MsBackendOfflineSql-method
#' @aliases filterMsLevel,MsBackendOfflineSql-method
#' @aliases filterPrecursorMzRange,MsBackendOfflineSql-method
#' @aliases filterPrecursorMzValues,MsBackendOfflineSql-method
#' @aliases filterRt,MsBackendOfflineSql-method
#' @aliases peaksData,MsBackendOfflineSql-method
#' @aliases peaksVariables,MsBackendOfflineSql-method
#' @aliases reset,MsBackendOfflineSql-method
#' @aliases show,MsBackendOfflineSql-method
#' @aliases spectraData,MsBackendOfflineSql-method
#' @aliases supportsSetBackend,MsBackendOfflineSql-method
#' @aliases tic,MsBackendOfflineSql-method
#' @aliases uniqueMsLevels,MsBackendOfflineSql-method
#' @aliases intensity,MsBackendOfflineSql-method
#' @aliases mz,MsBackendOfflineSql-method
#'
#' @description
#'
#' The `MsBackendSql` is an implementation for the [Spectra::MsBackend()] class
#' for [Spectra::Spectra()] objects which stores and retrieves MS data from a
#' SQL database. New databases can be created from raw MS data files using
#' `createMsBackendSqlDatabase()`.
#'
#' @details
#'
#' The `MsBackendSql` class is principally a *read-only* backend but by
#' extending the [Spectra::MsBackendCached()] backend from the `Spectra`
#' package it allows changing and adding (**temporarily**) spectra variables
#' **without** changing the original data in the SQL database.
#'
#' @note
#'
#' The `MsBackendSql` backend keeps an (open) connection to the SQL database
#' with the data and hence does not support saving/loading of a backend to
#' disk (e.g. using `save` or `saveRDS`). Also, for the same reason, the
#' `MsBackendSql` does not support parallel processing. The `backendBpparam()`
#' method for `MsBackendSql` will thus always return a
#' [BiocParallel::SerialParam()] object.
#'
#' The [MsBackendOfflineSql()] could be used as an alternative as it supports
#' saving/loading the data to/from disk and supports also parallel processing.
#'
#' @section Creation of backend objects:
#'
#' New backend objects can be created with the `MsBackendSql()` function.
#' SQL databases can be created and filled with MS data from raw data files
#' using the `createMsBackendSqlDatabase()` function or using
#' `backendInitialize()` and providing all data with parameter `data`. In
#' addition it is possible to create a database from a `Spectra` object
#' changing its backend to a `MsBackendSql` or `MsBackendOfflineSql` using
#' the [Spectra::setBackend()] function.
#' Existing SQL databases (created previously with
#' `createMsBackendSqlDatabase()` or `backendInitialize()` with the `data`
#' parameter) can be loaded using the *conventional* way to create/initialize
#' `MsBackend` classes, i.e. using `backendInitialize()`.
#'
#' - `createMsBackendSqlDatabase()`: create a database and fill it with MS data.
#'   Parameter `dbcon` is expected to be a database connection, parameter `x`
#'   a `character` vector with the file names from which to import the data.
#'   Parameter `backend` is used for the actual data import and defaults to
#'   `backend = MsBackendMzR()` hence allowing to import data from mzML, mzXML
#'   or netCDF files. Parameter `chunksize` allows to define the number of
#'   files (`x`) from which the data should be imported in one iteration. With
#'   the default `chunksize = 10L` data is imported from 10 files in `x` at
#'   the same time (if `backend` supports it even in parallel) and this data
#'   is then inserted into the database. Larger chunk sizes will require more
#'   memory and also larger disk space (as data import is performed through
#'   temporary files) but might eventually be faster. Parameter `blob` allows
#'   to define whether m/z and intensity values from a spectrum should be
#'   stored as a *BLOB* SQL data type in the database (`blob = TRUE`, the
#'   default) or if individual m/z and intensity values for each peak should
#'   be stored separately (`blob = FALSE`). The latter case results in a much
#'   larger database and slower performance of the `peaksData` function, but
#'   would allow to define custom (manual) SQL queries on individual peak
#'   values. For `blob = TRUE`, the peaks data can be stored in two different
#'   ways which can be selected with the additional parameter
#'   `peaksStorageMode`. The default `peaksStorageMode = "blob2"` stores the
#'   full peaks matrix of a spectrum as a single entry to the database (into a
#'   single database table column) while `peaksStorageMode = "blob"` (which
#'   was the default until version 1.7.2) stores the m/z and intensity
#'   vectors as BLOB data type into two separate database table columns.
#'   Performance for `peaksData()` is thus about twice as fast for
#'   `peaksStorageMode = "blob2"`.
#'   While data can be stored in any SQL database, at present it is suggested
#'   to use MySQL/MariaDB databases. For `dbcon` being a connection to a
#'   MySQL/MariaDB database, the tables will use the *ARIA* engine providing
#'   faster data access and will use *table partitioning*: tables are
#'   splitted into multiple partitions which can improve data insertion and
#'   index generation. Partitioning can be defined with the parameters
#'   `partitionBy` and `partitionNumber`. By default `partitionBy = "none"`
#'   no partitioning is performed. For `blob = TRUE` partitioning is usually
#'   not required. Only for `blob = FALSE ` and very large datasets it is
#'   suggested to enable table partitioning by selecting either
#'   `partitionBy = "spectrum"` or `partitionBy = "chunk"`. The first option
#'   assignes consecutive spectra to different partitions while the latter
#'   puts spectra from files part of the same *chunk* into the same partition.
#'   Both options have about the same performance but
#'   `partitionBy = "spectrum"` requires less disk space.
#'   Note that, while inserting the data takes a considerable amount of
#'   time, also the subsequent creation of database indices can take very
#'   long (even longer than data insertion for `blob = FALSE`).
#'
#' - `backendInitialize()`: get access and initialize a `MsBackendSql` object.
#'   Parameter `object` is supposed to be a `MsBackendSql` instance, created
#'   e.g. with `MsBackendSql()`. Parameter `dbcon` is expected to be a
#'   connection to an existing *MsBackendSql* SQL database (created e.g. with
#'   `createMsBackendSqlDatabase()`). To initialize a `MsBackendOfflineSql()`
#'   all information required for a [DBI::dbConnect()] call to connect to a
#'   database need to be provided.
#'   `backendInitialize()` can alternatively also be used to create a **new**
#'   `MsBackendSql` database using the optional
#'   `data` parameter. In this case, `dbcon` is expected to be a writeable
#'   connection to an empty database and `data` a `DataFrame` with the **full**
#'   spectra and peaks data to be inserted into this database. The format of
#'   `data` should match the format of the `DataFrame` returned by the
#'   `spectraData()` function and requires columns `"mz"` and `"intensity"`
#'   with the m/z and intensity values of each spectrum.
#'   The `backendInitialize()` call will then create all necessary tables in
#'   the database, will fill these tables with the provided data and will
#'   return an `MsBackendSql` for this database.
#'   The `MsBackendSql` and `MsBackendOfflineSql` objects
#'   support the [Spectra::setBackend()] method from `Spectra` to change
#'   from (any) backend to a `MsBackendSql`. Any parameters to the
#'   `backendInitialize()` function can be passed to `setBackend()`. Note
#'   however that chunk-wise (or parallel) processing needs to be disabled
#'   in this case by passing eventually `f = factor()` to the `setBackend()`
#'   call.
#'
#' - `supportsSetBackend()`: whether `MsBackendSql` supports the `setBackend()`
#'   method to change the `MsBackend` of a `Spectra` object to a
#'   `MsBackendSql`. Returns `TRUE`, thus, changing the backend to a
#'   `MsBackendSql` is supported **if** a writeable database connection
#'   is provided in addition with parameter `dbcon` (i.e.
#'   `setBackend(sps, MsBackendSql(), dbcon = con)` with `con` being a
#'   connection to an **empty** database would store the full spectra
#'   data from the `Spectra` object `sps` into the specified database and
#'   would return a `Spectra` object that uses a `MsBackendSql`).
#'
#' - `backendBpparam()`: whether a `MsBackendSql` supports parallel processing.
#'   Takes a `MsBackendSql` and a parallel processing setup (see
#'   [BiocParallel::bpparam()] for details) as input and always returns a
#'   [BiocParallel::SerialParam()] since
#'   `MsBackendSql` does **not** support parallel processing.
#'
#' - `dbconn()`: returns the connection to the database.
#'
#' @section Subsetting, merging and filtering data:
#'
#' `MsBackendSql` objects can be subsetted using the `[` or `extractByIndex()`
#' functions. Internally, this will simply subset the `integer` vector of the
#' primary keys and eventually cached data. The original data in the database
#' **is not** affected by any subsetting operation. Any subsetting operation
#' can be *undone* by resetting the object with the `reset()` function.
#' Subsetting in arbitrary order as well as index replication is supported.
#'
#' Multiple `MsBackendSql` objects can also be merged (combined) with the
#' `backendMerge()` function. Note that this requires that all `MsBackendSql`
#' objects are connected to the **same** database. This function is thus
#' mostly used for combining `MsBackendSql` objects that were previously
#' splitted using e.g. `split()`.
#'
#' In addition, `MsBackendSql` supports all other filtering methods available
#' through [Spectra::MsBackendCached()]. Implementation of filter functions
#' optimized for `MsBackendSql` objects are:
#'
#' - `filterDataOrigin()`: filter the object retaining spectra with `dataOrigin`
#'   spectra variable values matching the provided ones with parameter
#'   `dataOrigin`. The function returns the results in the order of the
#'   values provided with parameter `dataOrigin`.
#'
#' - `filterMsLevel()`: filter the object based on the MS levels specified with
#'   parameter `msLevel`. The function does the filtering using SQL queries.
#'   If `"msLevel"` is a *local* variable stored within the object (and hence
#'   in memory) the default implementation in `MsBackendCached` is used
#'   instead.
#'
#' - `filterPrecursorMzRange()`: filters the data keeping only spectra with a
#'   `precursorMz` within the m/z value range provided with parameter `mz`
#'   (i.e. all spectra with a precursor m/z `>= mz[1L]` and `<= mz[2L]`).
#'
#' - filterPrecursorMzValues()`: filters the data keeping only spectra with
#'   precursor m/z values matching the value(s) provided with parameter `mz`.
#'   Parameters `ppm` and `tolerance` allow to specify acceptable differences
#'   between compared values. Lengths of `ppm` and `tolerance` can be either
#'   `1` or equal to `length(mz)` to use different values for ppm and
#'   tolerance for each provided m/z value.
#'
#' - `filterRt()`: filter the object keeping only spectra with retention times
#'   within the specified retention time range (parameter `rt`). Optional
#'   parameter `msLevel.` allows to restrict the retention time filter only
#'   on the provided MS level(s) returning all spectra from other MS levels.
#'
#' @section Accessing and *modifying* data:
#'
#' The functions listed here are specifically implemented for `MsBackendSql`.
#' In addition, `MsBackendSql` inherits and supports all data accessor,
#' filtering functions and data manipulation functions from
#' [Spectra::MsBackendCached()].
#'
#' - `$`, `$<-`: access or set (add) spectra variables in `object`. Spectra
#'   variables added or modified using the `$<-` are *cached* locally within
#'   the object (data in the database is never changed). To restore an object
#'   (i.e. drop all cached values) the `reset` function can be used.
#'
#' - `dataStorage()`: returns a `character` vector same length as there are
#'   spectra in `object` with the name of the database containing the data.
#'
#' - `intensity<-`: not supported.
#'
#' - `mz<-`: not supported.
#'
#' - `peaksData()`: returns a `list` with the spectras' peak data. The length of
#'   the list is equal to the number of spectra in `object`. Each element of
#'   the list is a `matrix` with columns according to parameter `columns`. For
#'   an empty spectrum, a `matrix` with 0 rows is returned. Use
#'   `peaksVariables(object)` to list supported values for parameter
#'   `columns`.
#'
#' - `peaksVariables()`: returns a `character` with the available peak
#'   variables, i.e. columns that could be queried with `peaksData()`.
#'
#' - `reset()`: *restores* an `MsBackendSql` by re-initializing it with the
#'   data from the database. Any subsetting or cached spectra variables will
#'   be lost.
#'
#' - `spectraData()`: gets general spectrum metadata. `spectraData()` returns
#'   a `DataFrame` with the same number of rows as there are spectra in
#'   `object`. Parameter `columns` allows to select specific spectra
#'   variables.
#'
#' - `spectraNames()`, `spectraNames<-`: returns a `character` of length equal
#'   to the number of spectra in `object` with the primary keys of the spectra
#'   from the database (converted to `character`). Replacing spectra names
#'   with `spectraNames<-` is not supported.
#'
#' - `uniqueMsLevels()`: returns the unique MS levels of all spectra in
#'   `object`.
#'
#' - `tic()`: returns the originally reported total ion count (for
#'   `initial = TRUE`) or calculates the total ion count from the intensities
#'   of each spectrum (for `initial = FALSE`).
#'
#' @section Implementation notes:
#'
#' Internally, the `MsBackendSql` class contains only the primary keys for all
#' spectra stored in the SQL database. Keeping only these `integer` in memory
#' guarantees a minimal memory footpring of the object. Still, depending of
#' the number of spectra in the database, this `integer` vector might become
#' very large. Any data access will involve SQL calls to retrieve the data
#' from the database. By extending the [Spectra::MsBackendCached()] object
#' from the `Spectra` package, the `MsBackendSql` supports to (temporarily,
#' i.e. for the duration of the R session) add or modify spectra variables.
#' These are however stored in a `data.frame` within the object thus
#' increasing the memory demand of the object.
#'
#' @param backend For `createMsBackendSqlDatabase()`: MS backend that can be
#'     used to import MS data from the raw files specified with
#'     parameter `x`.
#'
#' @param blob For `createMsBackendSqlDatabase()`, `setBackend()`:
#'     `logical(1)` whether individual m/z and intensity values should be
#'     stored separately (`blob = FALSE`) or if the peaks data should be
#'     stored as data type *BLOB* into the database (`blob = TRUE`, the
#'     default). See also parameter `peaksStorageMode` for different data
#'     storage options.
#'
#' @param BPPARAM for `backendBpparam()`: `BiocParallel` parallel processing
#'     setup. See [BiocParallel::bpparam()] for more information.
#'
#' @param chunksize For `createMsBackendSqlDatabase()`: `integer(1)` defining
#'     the number of input that should be processed per iteration. With
#'     `chunksize = 1` each file specified with `x` will be imported and its
#'     data inserted to the database. With `chunksize = 5` data from 5 files
#'     will be imported (in parallel) and inserted to the database. Thus,
#'     higher values might result in faster database creation, but require
#'     also more memory.
#'
#' @param columns For `spectraData()`: `character()` optionally defining a
#'     subset of spectra variables that should be returned. Defaults to
#'     `columns = spectraVariables(object)` hence all variables are returned.
#'     For `peaksData` accessor: optional `character` with requested columns
#'     in the individual `matrix` of the returned `list`. Defaults to
#'     `columns = c("mz", "intensity")` but all columns listed by
#'     `peaksVariables` would be supported.
#'
#' @param data For `backendInitialize()`: optional `DataFrame` with the full
#'     spectra data that should be inserted into a (new) `MsBackendSql`
#'     database. If provided, it is assumed that `dbcon` is a (writeable)
#'     connection to an empty database into which `data` should be inserted.
#'     `data` could be the output of `spectraData` from another backend.
#'
#' @param dataOrigin For `filterDataOrigin()`: `character` with *data origin*
#'     values to which the data should be subsetted.
#'
#' @param dbcon Connection to a database.
#'
#' @param drop For `[`: `logical(1)`, ignored.
#'
#' @param initial For `tic()`: `logical(1)` whether the original total ion count
#'     should be returned (`initial = TRUE`, the default) or whether it
#'     should be calculated on the spectras' intensities (`initial = FALSE`).
#'
#' @param i For `[`: `integer` or `logical` to subset the object.
#'
#' @param j For `[`: ignored.
#'
#' @param msLevel For `filterMsLevel()`: `integer` specifying the MS levels to
#'     filter the data.
#'
#' @param msLevel. For `filterRt(): `integer` with the MS level(s) on which the
#'     retention time filter should be applied (all spectra from other MS
#'     levels are considered for the filter and are returned *as is*). If not
#'     specified, the retention time filter is applied to all MS levels in
#'     `object`.
#'
#' @param mz For `filterPrecursorMzRange()`: `numeric(2)` with the desired lower
#'     and upper limit of the precursor m/z range.
#'     For `filterPrecursorMzValues()`: `numeric` with the m/z value(s) to
#'     filter the object.
#'
#' @param name For `<-`: `character(1)` with the name of the spectra variable
#'     to replace.
#'
#' @param object A `MsBackendSql` instance.
#'
#' @param partitionBy For `createMsBackendSqlDatabase()`: `character(1)`
#'     defining if and how the peak data table should be partitioned. `"none"`
#'     (default): no partitioning, `"spectrum"`: peaks are assigned to the
#'     partition based on the spectrum ID (number), i.e. spectra are evenly
#'     (consecutively) assigned across partitions. For `partitionNumber = 3`,
#'     the first spectrum is assigned to the first partition, the second to
#'     the second, the third to the third and the fourth spectrum again to
#'     the first partition. `"chunk"`: spectra processed as part of the same
#'     *chunk* are placed into the same partition. All spectra from the next
#'     processed chunk are assigned to the next partition. Note that this is
#'     only available for MySQL/MariaDB databases, i.e., if `con` is a
#'     `MariaDBConnection`.
#'     See details for more information.
#'
#' @param partitionNumber For `createMsBackendSqlDatabase()`: `integer(1)`
#'     defining the number of partitions the database table will be
#'     partitioned into (only supported for MySQL/MariaDB databases).
#'
#' @param peaksStorageMode `character(1)` defining how peaks variables are
#'     stored in the database. The default `peaksStorageMode = "blob2"` stores
#'     the full peaks matrix of each spectrum as data type *BLOB* as one entry
#'     into a single database table column. `peaksStorageMode = "blob"` stores
#'     in contrast the m/z and intensity vectors as separate BLOB types into
#'     two database tables. `peaksStorageMode = "long"` allows to store the
#'     data in *long mode*, i.e. each intensity and m/z value is stored
#'     individually in the database.
#'
#' @param ppm For `filterPrecursorMzValues()`: `numeric` with the m/z-relative
#'     maximal acceptable difference for a m/z value to be considered
#'     matching. Can be of length 1 or equal to `length(mz)`.
#'
#' @param rt For `filterRt()`: `numeric(2)` with the lower and upper retention
#'     time. Spectra with a retention time `>= rt[1]` and `<= rt[2]` are
#'     returned.
#'
#' @param tolerance For `filterPrecursorMzValues()`: `numeric` with the absolute
#'     difference for m/z values to be considered matching. Can be of length 1
#'     or equal to `length(mz)`.
#'
#' @param value For all setter methods: replacement value.
#'
#' @param x For `createMsBackendSqlDatabase()`: `character` with the names of
#'     the raw data files from which the data should be imported. For other
#'     methods an `MsqlBackend` instance.
#'
#' @param ... For `[`: ignored. For `backendInitialize()`, if parameter `data`
#'     is used: additional parameters to be passed to the function creating the
#'     database such as `blob` or `peaksStorageMode`. For `setBackend()`: any
#'     parameters supported by `backendInitialize()` or
#'     `createMsBackendSqlDatabase()`.
#'
#' @name MsBackendSql
#'
#' @return See documentation of respective function.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @exportClass MsBackendSql
#'
#' @examples
#'
#' ####
#' ## Create a new MsBackendSql database
#'
#' ## Define a file from which to import the data
#' data_file <- system.file("microtofq", "MM8.mzML", package = "msdata")
#'
#' ## Create a database/connection to a database
#' library(RSQLite)
#' db_file <- tempfile()
#' dbc <- dbConnect(SQLite(), db_file)
#'
#' ## Import the data from the file into the database
#' createMsBackendSqlDatabase(dbc, data_file)
#' dbDisconnect(dbc)
#'
#' ## Initialize a MsBackendSql
#' dbc <- dbConnect(SQLite(), db_file)
#' be <- backendInitialize(MsBackendSql(), dbc)
#'
#' be
#'
#' ## Original data source
#' head(be$dataOrigin)
#'
#' ## Data storage
#' head(dataStorage(be))
#'
#' ## Access all spectra data
#' spd <- spectraData(be)
#' spd
#'
#' ## Available variables
#' spectraVariables(be)
#'
#' ## Access mz values
#' mz(be)
#'
#' ## Subset the object to spectra in arbitrary order
#' be_sub <- be[c(5, 1, 1, 2, 4, 100)]
#' be_sub
#'
#' ## The internal spectrum IDs (primary keys from the database)
#' be_sub$spectrum_id_
#'
#' ## Add additional spectra variables
#' be_sub$new_variable <- "B"
#'
#' ## This variable is *cached* locally within the object (not inserted into
#' ## the database)
#' be_sub$new_variable
NULL

#' @importClassesFrom DBI DBIConnection
#'
#' @noRd
setClassUnion("DBIConnectionOrNULL", c("DBIConnection", "NULL"))

#' @importClassesFrom S4Vectors DataFrame
#'
#' @importClassesFrom Spectra MsBackendCached
setClass(
    "MsBackendSql",
    contains = "MsBackendCached",
    slots = c(
        dbcon = "DBIConnectionOrNULL",
        spectraIds = "integer",
        .tables = "list",
        peak_fun = "function"),
    prototype = prototype(
        dbcon = NULL,
        spectraIds = integer(),
        .tables = list(),
        peak_fun = .fetch_peaks_data_long,
        readonly = TRUE, version = "0.2"))

#' @importFrom methods .valueClassTest is new validObject
#'
#' @noRd
setValidity("MsBackendSql", function(object) {
    msg <- NULL
    if (length(object@spectraIds) != object@nspectra)
        msg <- stri_c("Number of spectra IDs does not match the number ",
                      "of spectra")
    ## Can not validate the connection because of the MsBackendOfflineSql
    ## msg <- .valid_dbcon(object@dbcon)
    if (is.null(msg)) TRUE
    else msg
})

#' @importMethodsFrom Spectra show
#'
#' @exportMethod show
#'
#' @rdname MsBackendSql
setMethod("show", "MsBackendSql", function(object) {
    callNextMethod()
    if (!is.null(.dbcon(object))) {
        info <- dbGetInfo(.dbcon(object))
        cat("Database: ", info$dbname, "\n", sep = "")
    }
})

#' @exportMethod backendInitialize
#'
#' @importMethodsFrom ProtGenerics backendInitialize
#'
#' @importFrom DBI dbGetQuery
#'
#' @rdname MsBackendSql
setMethod("backendInitialize", "MsBackendSql",
          function(object, dbcon, data, ...) {
    if (missing(dbcon))
        stop("Parameter 'dbcon' is required for 'MsBackendSql'")
    ## if data is provided create a new database in dbcon with the data.
    if (!missing(data))
        .create_from_spectra_data(dbcon, data = data, ...)
    msg <- .valid_dbcon(dbcon)
    if (length(msg)) stop(msg)
    object@dbcon <- dbcon
    object@spectraIds <- dbGetQuery(
        dbcon, stri_c("select spectrum_id_ from msms_spectrum ",
                      "order by spectrum_id_"))[, 1L]
    object@.tables <- list(
        msms_spectrum = colnames(
            dbGetQuery(dbcon, "select * from msms_spectrum limit 0")))
    ## Whether  m/z and intensity values are stored as BLOBs
    if (any(dbListTables(dbcon) == "msms_spectrum_peak_blob"))
        object@peak_fun <- .fetch_peaks_data_blob
    if (any(dbListTables(dbcon) == "msms_spectrum_peak_blob2"))
        object@peak_fun <- .fetch_peaks_data_blob2
    ## Initialize cached backend
    object <- callNextMethod(
        object, nspectra = length(object@spectraIds),
        spectraVariables = c(unique(unlist(object@.tables))))
    validObject(object)
    object
})

#' @exportMethod dataStorage
#'
#' @importMethodsFrom ProtGenerics dataStorage
#'
#' @importFrom DBI dbGetInfo
#'
#' @rdname MsBackendSql
setMethod("dataStorage", "MsBackendSql", function(object) {
    if (!is.null(.dbcon(object))) {
        info <- dbGetInfo(.dbcon(object))
        rep(info$dbname, length(object))
    } else character()
})

#' @exportMethod [
#'
#' @importFrom MsCoreUtils i2index
#'
#' @importFrom methods slot<-
#'
#' @importFrom S4Vectors extractROWS
#'
#' @rdname MsBackendSql
setMethod("[", "MsBackendSql", function(x, i, j, ..., drop = FALSE) {
    if (missing(i))
        return(x)
    i <- i2index(i, length(x), x@spectraIds)
    extractByIndex(x, i)
})

#' @rdname MsBackendSql
#'
#' @importMethodsFrom ProtGenerics extractByIndex
#'
#' @export
setMethod("extractByIndex", c("MsBackendSql", "ANY"), function(object, i) {
    slot(object, "spectraIds", check = FALSE) <- object@spectraIds[i]
    callNextMethod(object, i = i)
})

#' @importMethodsFrom ProtGenerics peaksData
#'
#' @exportMethod peaksData
#'
#' @rdname MsBackendSql
setMethod(
    "peaksData", "MsBackendSql",
    function(object, columns = c("mz", "intensity")) {
        object@peak_fun(object, columns)
    })

#' @importMethodsFrom ProtGenerics peaksVariables
#'
#' @exportMethod peaksVariables
#'
#' @rdname MsBackendSql
setMethod("peaksVariables", "MsBackendSql",
          function(object) .available_peaks_variables(object))

#' @exportMethod intensity
#'
#' @importMethodsFrom ProtGenerics intensity
#'
#' @rdname MsBackendSql
setMethod("intensity", "MsBackendSql", function(object) {
    NumericList(object@peak_fun(object, "intensity", TRUE), compress = FALSE)
})

#' @exportMethod intensity<-
#'
#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname MsBackendSql
setReplaceMethod("intensity", "MsBackendSql", function(object, value) {
    stop("Can not replace original intensity values in the database.")
})

#' @exportMethod mz
#'
#' @importMethodsFrom ProtGenerics mz
#'
#' @rdname MsBackendSql
setMethod("mz", "MsBackendSql", function(object) {
    NumericList(object@peak_fun(object, "mz", TRUE), compress = FALSE)
})

#' @exportMethod mz<-
#'
#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname MsBackendSql
setReplaceMethod("mz", "MsBackendSql", function(object, value) {
    stop("Can not replace original data  in the database.")
})

#' @rdname MsBackendSql
#'
#' @export
setReplaceMethod("$", "MsBackendSql", function(x, name, value) {
    if (name %in% c("spectrum_id_"))
        stop("Spectra IDs can not be changed.", call. = FALSE)
    callNextMethod()
})

#' @importMethodsFrom ProtGenerics spectraData spectraVariables
#'
#' @exportMethod spectraData
#'
#' @rdname MsBackendSql
setMethod("spectraData", "MsBackendSql",
          function(object, columns = spectraVariables(object)) {
              .spectra_data_sql(object, columns = columns)
          })

#' @exportMethod reset
#'
#' @importMethodsFrom Spectra reset
#'
#' @rdname MsBackendSql
setMethod("reset", "MsBackendSql", function(object) {
    message("Restoring original data ...", appendLF = FALSE)
    if (is(object@dbcon, "DBIConnection"))
        object <- backendInitialize(MsBackendSql(), object@dbcon)
    message("DONE")
    object
})

#' @exportMethod spectraNames
#'
#' @importMethodsFrom ProtGenerics spectraNames
#'
#' @rdname MsBackendSql
setMethod("spectraNames", "MsBackendSql", function(object) {
    as.character(object@spectraIds)
})

#' @exportMethod spectraNames<-
#'
#' @importMethodsFrom ProtGenerics spectraNames<-
#'
#' @rdname MsBackendSql
setReplaceMethod("spectraNames", "MsBackendSql", function(object, value) {
    stop("Replacing spectraNames is not supported for ", class(object)[1L])
})

#' @importMethodsFrom ProtGenerics filterMsLevel uniqueMsLevels
#'
#' @rdname MsBackendSql
#'
#' @exportMethod filterMsLevel
setMethod(
    "filterMsLevel", "MsBackendSql",
    function(object, msLevel = uniqueMsLevels(object)) {
        if (!length(msLevel))
            return(object)
        if(.has_local_variable(object, "msLevel"))
            callNextMethod()
        else {
            qry <- stri_c(.id_query(object),
                          "msLevel in (", stri_c(msLevel, collapse = ","),")")
            .subset_query(object, qry)
        }
    })

#' @importMethodsFrom ProtGenerics filterRt
#'
#' @rdname MsBackendSql
#'
#' @exportMethod filterRt
setMethod("filterRt", "MsBackendSql", function(object, rt = numeric(),
                                               msLevel. = integer()) {
    if (!length(rt) || all(is.infinite(rt)))
        return(object)
    rt <- range(rt)
    if (.has_local_variable(object, "rtime") |
        (length(msLevel.) && .has_local_variable(object, "msLevel"))) {
        if (length(msLevel.)) {
            if (!.has_local_variable(object, "msLevel"))
                object$msLevel <- object$msLevel
        }
        if (!.has_local_variable(object, "rtime"))
            object$rtime <- object$rtime
        callNextMethod()
    } else {
        if (rt[1L] == -Inf)
            rt[1L] <- -1e12
        if (rt[2L] == Inf)
            rt[2L] <- 1e12
        if (length(msLevel.)) {
            msl <- stri_c(msLevel., collapse = ",")
            qry <- stri_c(.id_query(object),
                          "(rtime >= ", rt[1L], " and rtime <= ", rt[2L],
                          " and msLevel in (", msl, ")) or msLevel not in (",
                          msl, ")")
        } else {
            qry <- stri_c(.id_query(object),
                          "rtime >= ", rt[1L], " and rtime <= ", rt[2L], "")
        }
        .subset_query(object, qry)
    }
})

#' @importMethodsFrom ProtGenerics filterDataOrigin dataOrigin
#'
#' @rdname MsBackendSql
#'
#' @exportMethod filterDataOrigin
setMethod(
    "filterDataOrigin", "MsBackendSql", function(object,
                                                 dataOrigin = character()) {
        if (!length(dataOrigin))
            return(object)
        if (.has_local_variable(object, "dataOrigin"))
            callNextMethod()
        else {
            qry <- stri_c(.id_query(object), "dataOrigin in (",
                          stri_c("'", dataOrigin, "'", collapse = ","),")")
            object <- .subset_query(object, qry)
            ## Need to ensure the order is correct.
            if (length(dataOrigin) > 1L)
                object <- extractByIndex(
                    object, order(fmatch(dataOrigin(object), dataOrigin)))
            object
        }
    })

#' @importMethodsFrom ProtGenerics filterPrecursorMzRange
#'
#' @rdname MsBackendSql
#'
#' @exportMethod filterPrecursorMzRange
setMethod(
    "filterPrecursorMzRange", "MsBackendSql", function(object,
                                                       mz = numeric()) {
        if (length(mz)) {
            if (.has_local_variable(object, "precursorMz"))
                callNextMethod()
            else {
                mz <- range(mz)
                qry <- stri_c(.id_query(object), "precursorMz >= ", mz[1L],
                              " and precursorMz <= ", mz[2L])
                .subset_query(object, qry)
            }
        } else object
    })

#' @importMethodsFrom ProtGenerics filterPrecursorMzValues
#'
#' @rdname MsBackendSql
#'
#' @exportMethod filterPrecursorMzValues
setMethod(
    "filterPrecursorMzValues", "MsBackendSql",
    function(object, mz = numeric(), ppm = 20, tolerance = 0) {
        if (length(mz)) {
            if (.has_local_variable(object, "precursorMz"))
                callNextMethod()
            else {
                qry <- stri_c(.id_query(object),
                              .precursor_mz_query(mz, ppm, tolerance))
                object <- .subset_query(object, qry)
                object
            }
        } else object
    })

#' @rdname MsBackendSql
#'
#' @exportMethod uniqueMsLevels
setMethod("uniqueMsLevels", "MsBackendSql", function(object, ...) {
    if (!is.null(.dbcon(object))) {
        dbGetQuery(.dbcon(object),
                   "select distinct msLevel from msms_spectrum")[, "msLevel"]
    } else integer()
})

#' @rdname MsBackendSql
#'
#' @importMethodsFrom ProtGenerics backendMerge
#'
#' @exportMethod backendMerge
setMethod("backendMerge", "MsBackendSql", function(object, ...) {
    object <- unname(c(object, ...))
    not_empty <- lengths(object) > 0
    if (any(not_empty))
        res <- .combine(object[not_empty])
    else res <- object[[1L]]
    validObject(res)
    res
})

#' @rdname MsBackendSql
#'
#' @importMethodsFrom ProtGenerics precScanNum
#'
#' @exportMethod precScanNum
setMethod("precScanNum", "MsBackendSql", function(object) {
    spectraData(object, "precScanNum")[, 1L]
})

#' @rdname MsBackendSql
#'
#' @importMethodsFrom ProtGenerics centroided
#'
#' @exportMethod centroided
setMethod("centroided", "MsBackendSql", function(object) {
    as.logical(callNextMethod())
})

#' @rdname MsBackendSql
#'
#' @importMethodsFrom ProtGenerics smoothed
#'
#' @exportMethod smoothed
setMethod("smoothed", "MsBackendSql", function(object) {
    as.logical(callNextMethod())
})

#' @importMethodsFrom ProtGenerics tic
#'
#' @importFrom Spectra intensity
#'
#' @importFrom MsCoreUtils vapply1d
#'
#' @exportMethod tic
#'
#' @rdname MsBackendSql
setMethod("tic", "MsBackendSql", function(object, initial = TRUE) {
    if (initial)
        spectraData(object, "totIonCurrent")[, 1L]
    else vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @importMethodsFrom Spectra supportsSetBackend
#'
#' @exportMethod supportsSetBackend
#'
#' @rdname MsBackendSql
setMethod("supportsSetBackend", "MsBackendSql", function(object, ...) {
    TRUE
})

#' @importMethodsFrom Spectra backendBpparam
#'
#' @importFrom BiocParallel SerialParam bpparam
#'
#' @rdname MsBackendSql
setMethod(
    "backendBpparam", signature = "MsBackendSql",
    function(object, BPPARAM = bpparam()) {
        SerialParam()
    })

#' @importMethodsFrom BiocGenerics dbconn
#'
#' @rdname MsBackendSql
setMethod("dbconn", "MsBackendSql", .dbcon)

#' @importMethodsFrom ProtGenerics setBackend
#'
#' @importFrom Spectra processingChunkFactor
#'
#' @noRd
setMethod(
    "setBackend", c("Spectra", "MsBackendSql"),
    function(object, backend, f = processingChunkFactor(object), dbcon, ...,
             BPPARAM = SerialParam()) {
        backend_class <- class(object@backend)[1L]
        if (missing(dbcon))
            stop("Parameter 'dbcon' is required for 'MsBackendSql'")
        if (length(object)) {
            .set_backend_insert_data(object, f = f, con = dbcon, ...)
            object@backend <- backendInitialize(backend, dbcon = dbcon)
        } else object@backend <- backendInitialize(
                   backend, data = spectraData(object@backend),
                   dbcon = dbcon, ...)
        object@processing <- Spectra:::.logging(object@processing,
                                                "Switch backend from ",
                                                backend_class, " to ",
                                                class(object@backend))
        object
    })
