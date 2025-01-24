#' @include MsBackendSql-functions.R MsBackendSql.R

#' @title SQL-based MS backend without active database connection
#'
#' @description
#'
#' The `MsBackendOfflineSql` backend extends the [MsBackendSql()] backend
#' directly and inherits thus all of its functions as well as properties.
#' The only difference between the two backend is that `MsBackendSql` keeps
#' an active connection to the SQL database inside the object while the
#' `MsBackendOfflineSql` backends reconnects to the SQL database for each
#' query. While the performance of the latter is slightly lower (due to the
#' need to connect/disconnect to the database for each function call) it can
#' also be used in a parallel processing environment.
#'
#' @section Creation of backend objects:
#'
#' An empty instance of an `MsBackendOfflineSql` class can be created using the
#' `MsBackendOfflineSql()` function. An existing *MsBackendSql* SQL database
#' can be loaded with the `backendInitialize()` function. This function takes
#' parameters `drv`, `dbname`, `user`, `password`, `host` and `port`, all
#' parameters that are passed to the `dbConnect()` function to connect to
#' the (**existing**) SQL database.
#'
#' See [MsBackendSql()] for information on how to create a *MsBackend* SQL
#' database.
#'
#' @param object A `MsBackendOfflineSql` object.
#'
#' @param data For `backendInitialize()`: optional `DataFrame` with the full
#'     spectra data that should be inserted into a (new) `MsBackendSql`
#'     database. If provided, it is assumed that the provided database
#'     connection information if for a (writeable) empty database into which
#'     `data` should be inserted. `data` could be the output of `spectraData`
#'     from another backend.
#'
#' @param drv A *DBI* database driver object (such as `SQLite()` from the
#'     `RSQLite` package or `MariaDB()` from the `RMariaDB` package). See
#'     [DBI::dbConnect()] for more information.
#'
#' @param dbname `character(1)` with the name of the database. Passed directly
#'     to [DBI::dbConnect()].
#'
#' @param user `character(1)` with the user name for the database. Passed
#'     directly to [DBI::dbConnect()].
#'
#' @param password `character(1)` with the password for the database. Note
#'     that this password is stored (unencrypted) within the object. Passed
#'     directly to [DBI::dbConnect()].
#'
#' @param host `character(1)` with the host running the database. Passed
#'     directly to [DBI::dbConnect()].
#'
#' @param port `integer(1)` with the port number (optional). Passed directly to
#'     [DBI::dbConnect()].
#'
#' @param ... ignored.
#'
#' @aliases MsBackendOfflineSql-class
#'
#' @name MsBackendOfflineSql
#'
#' @exportClass MsBackendOfflineSql
#'
#' @author Johannes Rainer
NULL

#' @importClassesFrom DBI DBIDriver
#'
#' @noRd
setClassUnion("DBIDriverOrNULL", c("DBIDriver", "NULL"))

setClass(
    "MsBackendOfflineSql",
    contains = "MsBackendSql",
    slots = c(
        driver = "DBIDriverOrNULL",
        dbname = "character",
        user = "character",
        password = "character",
        host = "character",
        port = "integer",
        flags = "integer"
    ),
    prototype = prototype(
        driver = NULL,
        dbcon = NULL,
        dbname = character(),
        user = character(),
        password = character(),
        host = character(),
        port = NA_integer_,
        flags = 1L,
        spectraIds = integer(),
        .tables = list(),
        readonly = TRUE, version = "0.1"))

#' @rdname MsBackendOfflineSql
#'
#' @export MsBackendOfflineSql
MsBackendOfflineSql <- function() {
    new("MsBackendOfflineSql")
}

#' @importMethodsFrom DBI dbConnect
.db_connect <- function(x) {
    if (!is.null(x@driver))
        dbConnect(x@driver, dbname = x@dbname, user = x@user,
                  password = x@password, host = x@host, port = x@port,
                  flags = x@flags)
    else NULL
}

#' @importMethodsFrom DBI dbDisconnect dbIsValid
.db_disconnect <- function(x) {
    if (!is.null(x@dbcon) && dbIsValid(x@dbcon))
        dbDisconnect(x@dbcon)
}

setMethod("show", "MsBackendOfflineSql", function(object) {
    object@dbcon <- .db_connect(object)
    callNextMethod()
    .db_disconnect(object)
})

#' @rdname MsBackendOfflineSql
setMethod("backendInitialize", "MsBackendOfflineSql",
          function(object, drv = NULL, dbname = character(),
                   user = character(), password = character(),
                   host = character(), port = NA_integer_, data, ...) {
              if (is.null(drv) || !inherits(drv, "DBIDriver"))
                  stop("Parameter 'drv' must be specified and needs to be ",
                       "an instance of 'DBIDriver' such as returned e.g. ",
                       "by 'SQLite()'")
              if (!length(dbname))
                  stop("At least the database name has to be provided with ",
                       "parameter 'dbname'. Other possibly required parameters",
                       " (depending on the used database system) could be ",
                       "'user', 'password', 'host' and 'port'.")
              object@driver <- drv
              object@dbname <- dbname
              object@user <- user
              object@password <- password
              object@host <- host
              object@port <- as.integer(port)
              object@flags <- 1L        # that's SQLITE_RO
              if (missing(data))
                  dbc <- dbConnect(drv, dbname = dbname, user = user,
                                   password = password, host = host,
                                   port = port, flags = object@flags)
              else
                  dbc <- dbConnect(drv, dbname = dbname, user = user,
                                   password = password, host = host,
                                   port = port)
              object <- callNextMethod(object, dbcon = dbc, data = data, ...)
              .db_disconnect(object)
              object
})

setMethod("dataStorage", "MsBackendOfflineSql",
          function(object) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("intensity", "MsBackendOfflineSql", function(object) {
    object@dbcon <- .db_connect(object)
    on.exit(.db_disconnect(object))
    callNextMethod()
})

setMethod("mz", "MsBackendOfflineSql", function(object) {
    object@dbcon <- .db_connect(object)
    on.exit(.db_disconnect(object))
    callNextMethod()
})

setMethod("peaksData", "MsBackendOfflineSql",
          function(object, columns = c("mz", "intensity")) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("peaksVariables", "MsBackendOfflineSql",
          function(object) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("spectraData", "MsBackendOfflineSql",
          function(object, columns = spectraVariables(object)) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("reset", "MsBackendOfflineSql",
          function(object) {
              message("Restoring original data ...", appendLF = FALSE)
              object <- backendInitialize(
                  MsBackendOfflineSql(), object@driver, dbname = object@dbname,
                  user = object@user, password = object@password,
                  host = object@host, port = object@port)
              message("DONE")
              object
          })

setMethod("filterMsLevel", "MsBackendOfflineSql",
          function(object, msLevel = uniqueMsLevels(object)) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("filterRt", "MsBackendOfflineSql",
          function(object, rt = numeric(), msLevel. = integer()) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("filterDataOrigin", "MsBackendOfflineSql",
          function(object, dataOrigin = character()) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("filterPrecursorMzRange", "MsBackendOfflineSql",
          function(object, mz = numeric()) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("filterPrecursorMzValues", "MsBackendOfflineSql",
          function(object, mz = numeric(), ppm = 20, tolerance = 0) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("uniqueMsLevels", "MsBackendOfflineSql",
          function(object, ...) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("backendMerge", "MsBackendOfflineSql",
          function(object, ...) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("tic", "MsBackendOfflineSql",
          function(object, initial = TRUE) {
              object@dbcon <- .db_connect(object)
              on.exit(.db_disconnect(object))
              callNextMethod()
          })

setMethod("backendBpparam", signature = "MsBackendOfflineSql",
          function(object, BPPARAM = bpparam()) {
              BPPARAM
          })

setMethod("dbconn", "MsBackendOfflineSql", .db_connect)

setMethod(
    "setBackend", c("Spectra", "MsBackendOfflineSql"),
    function(object, backend, f = processingChunkFactor(object), drv = NULL,
             dbname = character(), user = character(), password = character(),
             host = character(), port = NA_integer_, ...,
             BPPARAM = SerialParam()) {
        backend_class <- class(object@backend)[1L]
        if (is.null(drv))
            stop("Parameter 'drv' must be specified and needs to be ",
                 "an instance of 'DBIDriver' such as returned e.g. ",
                 "by 'SQLite()'")
        if (!length(dbname))
            stop("At least the database name has to be provided with ",
                 "parameter 'dbname'. Other possibly required parameters",
                 " (depending on the used database system) could be ",
                 "'user', 'password', 'host' and 'port'.")
        if (length(object)) {
            dbcon <- dbConnect(drv, dbname = dbname, user = user,
                               password = password, host = host,
                               port = port)
            .set_backend_insert_data(object, f = f, con = dbcon, ...)
            object@backend <- backendInitialize(
                backend, drv = drv, dbname = dbname, user = user,
                password = password, host = host, port = port)
            dbDisconnect(dbcon)
        } else object@backend <- backendInitialize(
                   backend, data = spectraData(object@backend), drv = drv,
                   dbname = dbname, user = user, password = password,
                   host = host, port = port, ...)
        object@processing <- Spectra:::.logging(object@processing,
                                                "Switch backend from ",
                                                backend_class, " to ",
                                                class(object@backend))
        object
})

## setReplaceMethod("$", "MsBackendOfflineSql", function(x, name, value) {
##     object@dbcon <- .db_connect(object)
##     on.exit(.db_disconnect(object))
##     res <- callNextMethod()
##     res
## })
