#' @title `Spectra` MS backend storing data in a SQL database
#'
#' @aliases MsqlBackend-class
#'
#' @description
#'
#' The `MsqlBackend` is an implementation for the [MsBackend()] class for
#' [Spectra()] objects which stores and retrieves MS data from a SQL database.
#' New databases can be created from raw MS data files using
#' `createMsqlBackendDatabase`.
#'
#' @details
#'
#' The `MsqlBackend` class is in principle a *read-only* backend but by
#' extending the [MsBackendCached()] backend from the `Spectra` package it
#' allows changing and adding (**temporarily**) spectra variables **without**
#' changing the original data in the SQL database.
#'
#' @section Creation of backend objects:
#'
#' SQL databases can be created and filled with MS data from raw data files
#' using the `createMsqlBackendDatabase` function. Existing SQL databases
#' (created previously with `createMsqlBackendDatabase` can be loaded using
#' the conventional way to create/initialize `MsBackend` classes, i.e. using
#' `backendInitialize`.
#'
#' - `createMsqlBackendDatabase`: create a database and fill it with MS data.
#'   Parameter `dbcon` is expected to be a database connection, parameter `x` a
#'   `character` vector with the file names from which to import the data.
#'   Parameter `backend` is used for the actual data import and defaults to
#'   `backend = MsBackendMzR()` hence allowing to import data from mzML, mzXML
#'   or netCDF files. Parameter `chunksize` allows to define the number of
#'   files (`x`) from which the data should be imported in one iteration. With
#'   the default `chunksize = 10L` data is imported from 10 files in `x` at the
#'   same time (if `backend` supports it even in parallel) and this data is then
#'   inserted into the database. Larger chunk sizes will require more memory and
#'   also larger disk space (as data import is performed through temporary
#'   files) but might eventually be faster. While data can be stored in any SQL
#'   database, at present it is suggested to use MySQL/MariaDB databases. For
#'   `dbcon` being a connection to a MySQL/MariaDB database, the tables will use
#'   the *ARIA* engine providing faster data access and will use table
#'   partitioning (using by default 10 partitions). Note that, while inserting
#'   the data takes a considerable amount of time, also the subsequent creation
#'   of database indices can take very long (even longer than data insertion).
#'
#' - `backendInitialize`: get access and initialize a `MsqlBackend` object.
#'   Parameter `object` is supposed to be a `MsqlBackend` instance, created e.g.
#'   with `MsqlBackend()`. Parameter `dbcon` is expected to be a connection to
#'   a SQL database previously created with the `createMsqlBackendDatabase`
#'   function.
#'
#' @section Subsetting and filtering data:
#'
#' `MsqlBackend` objects can be subsetted using the `[` function. Internally,
#' this will simply subset the `integer` vector of the primary keys and
#' eventually cached data. The original data in the database **is not** affected
#' by any subsetting operation. Any subsetting operation can be *undone* by
#' resetting the object with the `reset` function. Subsetting in arbitrary
#' order as well as index replication is supported.
#'
#' In addition, `MsqlBackend` supports all other filtering methods available
#' through [MsBackendCached()].
#'
#' @section Accessing and *modifying* data:
#'
#' The functions listed here are specifically implemented for `MsqlBackend`. In
#' addition, `MsqlBackend` inherits and supports all data accessor, filtering
#' functions and data manipulation functions from [MsBackendCached()].
#'
#' - `$`, `$<-`: access or set (add) spectra variables in `object`. Spectra
#'   variables added or modified using the `$<-` are *cached* locally within
#'   the object (data in the database is never changed). To restore an object
#'   (i.e. drop all cached values) the `reset` function can be used.
#'
#' - `dataStorage`: returns a `character` vector same length as there are
#'   spectra in `object` with the name of the database containing the data.
#'
#' - `intensity<-`: not supported.
#'
#' - `mz<-`: not supported.
#'
#' - `peaksData`: returns a `list` with the spectras' peak data. The length of
#'   the list is equal to the number of spectra in `object`. Each element of
#'   the list is a `matrix` with columns according to parameter `columns`. For
#'   an empty spectrum, a `matrix` with 0 rows is returned. Use
#'   `peaksVariables(object)` to list supported values for parameter `columns`.
#'
#' - `peaksVariables`: returns a `character` with the available peak variables,
#'   i.e. columns that could be queried with `peaksData`.
#'
#' - `reset`: *restores* an `MsqlBackend` by re-initializing it with the data
#'   from the database. Any subsetting or cached spectra variables will be lost.
#'
#' - `spectraData`: gets or general spectrum metadata.  `spectraData` returns
#'   a `DataFrame` with the same number of rows as there are spectra in
#'   `object`. Parameter `columns` allows to select specific spectra variables.
#'
#' - `spectraNames`, `spectraNames<-`: returns a `character` of length equal to
#'   the number of spectra in `object` with the primary keys of the spectra from
#'   the database (converted to `character`). Replacing spectra names with
#'   `spectraNames<-` is not supported.
#'
#' @section Implementation notes:
#'
#' Internally, the `MsqlBackend` class contains only the primary keys for all
#' spectra stored in the SQL database. Keeping only these `integer` in memory
#' guarantees a minimal memory footpring of the object. Still, depending of the
#' number of spectra in the database, this `integer` vector might become very
#' large. Any data access will involve SQL calls to retrieve the data from the
#' database. By extending the [MsBackendCached()] object from the `Spectra`
#' package, the `MsqlBackend` supports to (temporarily, i.e. for the duration
#' of the R session) add or modify spectra variables. These are however stored
#' in a `data.frame` within the object thus increasing the memory demand of the
#' object.
#'
#' @param dbcon Connection to a database.
#'
#' @param backend For `createMsqlBackendDatabase`: MS backend that can be used
#'     to import MS data from the raw files specified with parameter `x`.
#'
#' @param chunksize For `createMsqlBackendDatabase`: `integer(1)` defining the
#'     number of input that should be processed per iteration. With
#'     `chunksize = 1` each file specified with `x` will be imported and its
#'     data inserted to the database. With `chunksize = 5` data from 5 files
#'     will be imported (in parallel) and inserted to the database. Thus, higher
#'     values might result in faster database creation, but require also more
#'     memory.
#'
#' @param columns For `spectraData`: `character()` optionally defining a subset
#'     of spectra variables that should be returned. Defaults to
#'     `columns = spectraVariables(object)` hence all variables are returned.
#'     For `peaksData` accessor: optional `character` with requested columns in
#'     the individual `matrix` of the returned `list`. Defaults to
#'     `columns = c("mz", "intensity")` but all columns listed by
#'     `peaksVariables` would be supported.
#'
#' @param drop For `[`: `logical(1)`, ignored.
#'
#' @param i For `[`: `integer` or `logical` to subset the object.
#'
#' @param j For `[`: ignored.
#'
#' @param name For `<-`: `character(1)` with the name of the spectra variable
#'     to replace.
#'
#' @param object A `MsqlBackend` instance.
#'
#' @param value For all setter methods: replacement value.
#'
#' @param x For `createMsqlBackendDatabase`: `character` with the names of the
#'     raw data files from which the data should be imported. For other methods
#'     an `MsqlBaackend` instance.
#'
#' @param ... For `[`: ignored.
#'
#' @name MsqlBackend
#'
#' @return See documentation of respective function.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @exportClass MsqlBackend
#'
#' @examples
#'
#' ####
#' ## Create a new MsqlBackend database
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
#' createMsqlBackendDatabase(dbc, data_file)
#' dbDisconnect(dbc)
#'
#' ## Initialize a MsqlBackend
#' dbc <- dbConnect(SQLite(), db_file)
#' be <- backendInitialize(MsqlBackend(), dbc)
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
#' ## This variable is *cached* locally within the object (not inserted into the
#' ## database)
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
    "MsqlBackend",
    contains = "MsBackendCached",
    slots = c(
        dbcon = "DBIConnectionOrNULL",
        spectraIds = "integer",
        .tables = "list"),
    prototype = prototype(
        dbcon = NULL,
        spectraIds = integer(),
        .tables = list(),
        readonly = TRUE, version = "0.1"))

#' @importFrom methods .valueClassTest is new validObject
#'
#' @noRd
setValidity("MsqlBackend", function(object) {
    msg <- .valid_dbcon(object@dbcon)
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname MsqlBackend
#'
#' @importMethodsFrom Spectra show
#'
#' @exportMethod show
setMethod("show", "MsBackendCached", function(object) {
    callNextMethod()
    if (!is.null(.dbcon(object))) {
        info <- dbGetInfo(.dbcon(object))
        cat("Database: ", info$dbname, "\n", sep = "")
    }
})

#' @exportMethod backendInitialize
#'
#' @importMethodsFrom Spectra backendInitialize
#'
#' @importFrom DBI dbGetQuery
#'
#' @rdname MsqlBackend
setMethod("backendInitialize", "MsqlBackend",
          function(object, dbcon, ...) {
    if (missing(dbcon))
        stop("Parameter 'dbcon' is required for 'MsqlBackend'")
    msg <- .valid_dbcon(dbcon)
    if (length(msg)) stop(msg)
    object@dbcon <- dbcon

    res <- dbGetQuery(
        dbcon, "select spectrum_id_ from msms_spectrum")
    object@spectraIds <- res[, "spectrum_id_"]
    object@.tables <- list(
        msms_spectrum = colnames(
            dbGetQuery(dbcon, "select * from msms_spectrum limit 0")))
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
#' @rdname MsqlBackend
setMethod("dataStorage", "MsqlBackend", function(object) {
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
#' @rdname MsqlBackend
setMethod("[", "MsqlBackend", function(x, i, j, ..., drop = FALSE) {
    if (missing(i))
        return(x)
    i <- i2index(i, length(x), x@spectraIds)
    slot(x, "spectraIds", check = FALSE) <- x@spectraIds[i]
    x <- callNextMethod(x, i = i)
    x
})

#' @importMethodsFrom Spectra peaksData
#'
#' @exportMethod peaksData
#'
#' @rdname MsqlBackend
setMethod(
    "peaksData", "MsqlBackend",
    function(object, columns = c("mz", "intensity")) {
        pks <- .fetch_peaks_sql(object, columns)
        f <- as.factor(pks$spectrum_id_)        # using levels does not work because we can have duplicates
        pks <- unname(split.data.frame(pks, f)[as.character(object@spectraIds)])
        idx <- seq_along(columns) + 1
        lapply(pks, function(z) {
            if (nrow(z))
                as.matrix(z[, idx, drop = FALSE], rownames.force = FALSE)
            else matrix(NA_real_, ncol = length(columns), nrow = 0,
                        dimnames = list(character(), columns))
        })
    })

#' @importMethodsFrom Spectra peaksVariables
#'
#' @exportMethod peaksVariables
#'
#' @rdname MsqlBackend
setMethod("peaksVariables", "MsqlBackend",
          function(object) .available_peaks_variables(object))

#' @exportMethod intensity<-
#'
#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname MsqlBackend
setReplaceMethod("intensity", "MsqlBackend", function(object, value) {
    stop("Can not replace original intensity values in the database.")
})

#' @exportMethod mz<-
#'
#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname MsqlBackend
setReplaceMethod("mz", "MsqlBackend", function(object, value) {
    stop("Can not replace original data  in the database.")
})

#' @rdname MsqlBackend
#'
#' @export
setReplaceMethod("$", "MsqlBackend", function(x, name, value) {
    if (name %in% c("spectrum_id_"))
        stop("Spectra IDs can not be changed.", call. = FALSE)
    callNextMethod()
})

#' @importMethodsFrom Spectra spectraData spectraVariables
#'
#' @exportMethod spectraData
#'
#' @rdname MsqlBackend
setMethod("spectraData", "MsqlBackend",
          function(object, columns = spectraVariables(object)) {
              .spectra_data_sql(object, columns = columns)
          })

#' @exportMethod reset
#'
#' @importMethodsFrom Spectra reset
#'
#' @rdname MsqlBackend
setMethod("reset", "MsqlBackend", function(object) {
    message("Restoring original data ...", appendLF = FALSE)
    if (is(object@dbcon, "DBIConnection"))
        object <- backendInitialize(MsqlBackend(), object@dbcon)
    message("DONE")
    object
})

#' @exportMethod spectraNames
#'
#' @importMethodsFrom ProtGenerics spectraNames
#'
#' @rdname MsqlBackend
setMethod("spectraNames", "MsqlBackend", function(object) {
    as.character(object@spectraIds)
})

#' @exportMethod spectraNames<-
#'
#' @importMethodsFrom ProtGenerics spectraNames<-
#'
#' @rdname MsqlBackend
setReplaceMethod("spectraNames", "MsqlBackend",
                 function(object, value) {
                     stop(class(object)[1],
                          " does not support replacing spectra names (IDs).")
})
