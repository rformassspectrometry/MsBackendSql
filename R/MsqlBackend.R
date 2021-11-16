#' @title MS backend accessing the MassBank MySQL database
#'
#' @aliases MsBackendMassbankSql-class compounds
#'
#' @description
#'
#' The `MsBackendMassbankSql` provides access to mass spectrometry data from
#' [MassBank](https://massbank.eu/MassBank/) by directly accessing its
#' MySQL/MariaDb database. In addition it supports adding new spectra variables
#' or *locally* changing spectra variables provided by MassBank (without
#' changing the original values in the database).
#'
#' Note that `MsBackendMassbankSql` requires a local installation of the
#' MassBank database since direct database access is not supported for the
#' *main* MassBank instance.
#'
#' Also, some of the fields in the MassBank database are not directly compatible
#' with `Spectra`, such as the *collision energy* which is not available as a
#' numeric value. The collision energy as available in MassBank is reported as
#' spectra variable `"collision_energy_text"`. Also, precursor m/z values
#' reported for some spectra can not be converted to a `numeric` and hence `NA`
#' is reported with the spectra variable `precursorMz` for these spectra. The
#' variable `"precursor_mz_text"` can be used to get the *original* precursor
#' m/z reported in MassBank.
#'
#' @param dbcon For `backendInitialize,MsBackendMassbankSql`: SQL database
#'     connection to the MassBank (MariaDb) database.
#'
#' @param columns For `spectraData` accessor: optional `character` with column
#'     names (spectra variables) that should be included in the
#'     returned `DataFrame`. By default, all columns are returned.
#'
#' @param drop For `[`: not considered.
#'
#' @param initial For `tic`: `logical(1)` whether the initially
#'     reported total ion current should be reported, or whether the
#'     total ion current should be (re)calculated on the actual data
#'     (`initial = FALSE`).
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[`: not supported.
#'
#' @param name For `$` and `$<-`: the name of the spectra variable to return
#'     or set.
#'
#' @param object Object extending `MsBackendMassbankSql`.
#'
#' @param spectraVariables For `selectSpectraVariables`: `character` with the
#'     names of the spectra variables to which the backend should be subsetted.
#'
#' @param use.names For `lengths`: whether spectrum names should be used.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param x Object extending `MsBackendMassbankSql`.
#'
#' @param ... Additional arguments.
#'
#'
#' @section Supported Backend functions:
#'
#' The following functions are supported by the `MsBackendMassbankSqlMassbankDb`.
#'
#' - `[`: subset the backend. Only subsetting by element (*row*/`i`) is
#'   allowed
#'
#' - `$`, `$<-`: access or set/add a single spectrum variable (column) in the
#'   backend.
#'
#' - `acquisitionNum`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `peaksData` returns a `list` with the spectras' peak data. The length of
#'   the list is equal to the number of spectra in `object`. Each element of
#'   the list is a `matrix` with columns `"mz"` and `"intensity"`. For an empty
#'   spectrum, a `matrix` with 0 rows and two columns (named `mz` and
#'   `intensity`) is returned.
#'
#' - `backendInitialize`: initialises the backend by retrieving the IDs of all
#'   spectra in the database. Parameter `dbcon` with the connection to the
#'   MassBank MySQL database is required.
#'
#' - `dataOrigin`: gets a `character` of length equal to the number of spectra
#'   in `object` with the *data origin* of each spectrum. This could e.g. be
#'   the mzML file from which the data was read.
#'
#' - `dataStorage`: returns `"<MassBank>"` for all spectra.
#'
#' - `centroided`, `centroided<-`: gets or sets the centroiding
#'   information of the spectra. `centroided` returns a `logical`
#'   vector of length equal to the number of spectra with `TRUE` if a
#'   spectrum is centroided, `FALSE` if it is in profile mode and `NA`
#'   if it is undefined. See also `isCentroided` for estimating from
#'   the spectrum data whether the spectrum is centroided.  `value`
#'   for `centroided<-` is either a single `logical` or a `logical` of
#'   length equal to the number of spectra in `object`.
#'
#' - `collisionEnergy`, `collisionEnergy<-`: gets or sets the
#'   collision energy for all spectra in `object`. `collisionEnergy`
#'   returns a `numeric` with length equal to the number of spectra
#'   (`NA_real_` if not present/defined), `collisionEnergy<-` takes a
#'   `numeric` of length equal to the number of spectra in `object`. Note that
#'   the collision energy description from MassBank are provided as spectra
#'   variable `"collisionEnergyText"`.
#'
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the `list` is equal to the number of
#'   `spectra` in `object`.
#'
#' - `ionCount`: returns a `numeric` with the sum of intensities for
#'   each spectrum. If the spectrum is empty (see `isEmpty`),
#'   `NA_real_` is returned.
#'
#' - `isCentroided`: a heuristic approach assessing if the spectra in
#'   `object` are in profile or centroided mode. The function takes
#'   the `qtl` th quantile top peaks, then calculates the difference
#'   between adjacent m/z value and returns `TRUE` if the first
#'   quartile is greater than `k`. (See `Spectra:::.isCentroided` for
#'   the code.)
#'
#' - `isEmpty`: checks whether a spectrum in `object` is empty
#'   (i.e. does not contain any peaks). Returns a `logical` vector of
#'   length equal number of spectra.
#'
#' - `isolationWindowLowerMz`, `isolationWindowLowerMz<-`: gets or sets the
#'   lower m/z boundary of the isolation window.
#'
#' - `isolationWindowTargetMz`, `isolationWindowTargetMz<-`: gets or sets the
#'   target m/z of the isolation window.
#'
#' - `isolationWindowUpperMz`, `isolationWindowUpperMz<-`: gets or sets the
#'   upper m/z boundary of the isolation window.
#'
#' - `isReadOnly`: returns a `logical(1)` whether the backend is *read
#'   only* or does allow also to write/update data.
#'
#' - `length`: returns the number of spectra in the object.
#'
#' - `lengths`: gets the number of peaks (m/z-intensity values) per
#'   spectrum.  Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `0` is returned.
#'
#' - `msLevel`: gets the spectra's MS level. Returns an `integer`
#'   vector (of length equal to the number of spectra) with the MS
#'   level for each spectrum (or `NA_integer_` if not available).
#'
#' - `mz`: gets the mass-to-charge ratios (m/z) from the
#'   spectra. Returns a [NumericList()] or length equal to the number of
#'   spectra, each element a `numeric` vector with the m/z values of
#'   one spectrum.
#'
#' - `polarity`, `polarity<-`: gets or sets the polarity for each
#'   spectrum.  `polarity` returns an `integer` vector (length equal
#'   to the number of spectra), with `0` and `1` representing negative
#'   and positive polarities, respectively. `polarity<-` expects an
#'   integer vector of length 1 or equal to the number of spectra.
#'
#' - `precursorCharge`, `precursorIntensity`, `precursorMz`,
#'   `precScanNum`, `precAcquisitionNum`: get the charge (`integer`),
#'   intensity (`numeric`), m/z (`numeric`), scan index (`integer`)
#'   and acquisition number (`interger`) of the precursor for MS level
#'   2 and above spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `reset`: restores the backend to its original state, i.e. deletes all
#'   locally modified data and reinitializes the backend to the full data
#'   available in the database.
#'
#' - `rtime`, `rtime<-`: gets or sets the retention times for each
#'   spectrum (in seconds). `rtime` returns a `numeric` vector (length equal to
#'   the number of spectra) with the retention time for each spectrum.
#'   `rtime<-` expects a numeric vector with length equal to the
#'   number of spectra.
#'
#' - `scanIndex`: returns an `integer` vector with the *scan index*
#'   for each spectrum. This represents the relative index of the
#'   spectrum within each file. Note that this can be different to the
#'   `acquisitionNum` of the spectrum which is the index of the
#'   spectrum as reported in the mzML file.
#'
#' - `selectSpectraVariables`: reduces the information within the backend to
#'   the selected spectra variables.
#'
#' - `smoothed`,`smoothed<-`: gets or sets whether a spectrum is
#'   *smoothed*. `smoothed` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraData`: gets general spectrummetadata (annotation, also called
#'   header).  `spectraData` returns a `DataFrame`. Note that replacing the
#'   spectra data with `spectraData<-` is not supported.
#'
#' - `spectraNames`: returns a `character` vector with the names of
#'   the spectra in `object`.
#'
#' - `spectraVariables`: returns a `character` vector with the
#'   available spectra variables (columns, fields or attributes)
#'   available in `object`. This should return **all** spectra variables which
#'   are present in `object`, also `"mz"` and `"intensity"` (which are by
#'   default not returned by the `spectraVariables,Spectra` method).
#'
#' - `tic`: gets the total ion current/count (sum of signal of a
#'   spectrum) for all spectra in `object`. By default, the value
#'   reported in the original raw data file is returned. For an empty
#'   spectrum, `NA_real_` is returned.
#'
#' @section Not supported Backend functions:
#'
#' The following functions are not supported by the `MsBackendMassbankSql` since
#' the original data can not be changed.
#'
#' `backendMerge`, `export`, `filterDataStorage`, `filterPrecursorScan`,
#' `peaksData<-`, `filterAcquisitionNum`, `intensity<-`, `mz<-`, `precScanNum`,
#' `spectraData<-`, `spectraNames<-`.
#'
#' @section Retrieving compound annotations for spectra:
#'
#' While compound annotations are also provided *via* the `spectraVariables` of
#' the backend, it would also be possible to use the `compounds` function on
#' a `Spectra` object (that uses a `MsBackendMassbankSql` backend) to retrieve
#' compound annotations for the specific spectra.
#'
#' @name MsBackendMassbankSql
#'
#' @return See documentation of respective function.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @exportClass MsBackendMassbankSql
#'
#' @examples
#'
#' ## Create a connection to a database with MassBank data - in the present
#' ## example we connect to a tiny SQLite database bundled in this package
#' ## as public access to the MassBank MySQL is not (yet) supported. See the
#' ## vignette for more information on how to install MassBank locally and
#' ## enable MySQL database connections
#' library(RSQLite)
#' con <- dbConnect(SQLite(), system.file("sql", "minimassbank.sqlite",
#'     package = "MsBackendMassbank"))
#'
#' ## Given that we have the connection to a MassBank databas we can
#' ## initialize the backend:
#' be <- backendInitialize(MsBackendMassbankSql(), dbcon = con)
#' be
#'
#' ## Access MS level
#' msLevel(be)
#' be$msLevel
#'
#' ## Access m/z values
#' be$mz
#'
#' ## Access the full spectra data (including m/z and intensity values)
#' spectraData(be)
#'
#' ## Add a new spectra variable
#' be$new_variable <- "b"
#' be$new_variable
#'
#' ## Subset the backend
#' be_sub <- be[c(3, 1)]
#'
#' spectraNames(be)
#' spectraNames(be_sub)
NULL

setClassUnion("DBIConnectionOrNULL", c("DBIConnection", "NULL"))

#' @importClassesFrom DBI DBIConnection
#'
#' @importClassesFrom S4Vectors DataFrame
setClass(
    "MsqlBackend",
    contains = "MsBackend",
    slots = c(
        dbcon = "DBIConnectionOrNULL",
        spectraIds = "character",
        spectraVariables = "character",
        coreSpectraVariables = "character",
        localData = "DataFrame",
        .tables = "list"),
    prototype = prototype(
        dbcon = NULL,
        spectraIds = character(),
        spectraVariables = character(),
        coreSpectraVariables = names(Spectra:::.SPECTRA_DATA_COLUMNS),
        localData = DataFrame(),
        .tables = list(),
        readonly = TRUE, version = "0.1"))

#' @importFrom methods .valueClassTest is new validObject
#'
#' @noRd
setValidity("MsqlBackend", function(object) {
    msg <- .valid_dbcon(object@dbcon)
    msg <- c(msg, .valid_local_data(object@localData, object@spectraIds))
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname MsBackendMassbankSql
#'
#' @importFrom utils capture.output head tail
#'
#' @importMethodsFrom methods show
#'
#' @export
setMethod("show", "MsBackendMassbankSql", function(object) {
    n <- length(object@spectraIds)
    cat(class(object), "with", n, "spectra\n")
    if (n) {
        idx <- union(1:min(6, n), max(1, n-5):n)
        spd <- spectraData(object[idx, ],
                           c("msLevel", "precursorMz", "polarity"))
        if (!length(rownames(spd)))
            rownames(spd) <- idx
        txt <- capture.output(print(spd))
        cat(txt[-1], sep = "\n")
        sp_cols <- spectraVariables(object)
        cat(" ...", length(sp_cols) - 3, "more variables/columns.\n", "Use ",
            "'spectraVariables' to list all of them.\n")
    }
})

#' @exportMethod backendInitialize
#'
#' @importFrom DBI dbGetQuery
#' @importFrom S4Vectors make_zero_col_DFrame
#'
#' @rdname MsBackendMassbankSql
setMethod("backendInitialize", "MsBackendMassbankSql",
          function(object, dbcon, ...) {
    if (missing(dbcon))
        stop("Parameter 'dbcon' is required for 'MsBackendMassbankSql'")
    msg <- .valid_dbcon(dbcon)
    object@dbcon <- dbcon
    if (length(msg))
        stop(msg)

    res <- dbGetQuery(
        dbcon, "select spectrum_id, precursor_mz_text from msms_spectrum")
    object@spectraIds <- as.character(res[, "spectrum_id"])
    object@localData <- make_zero_col_DFrame(length(object@spectraIds))

    suppressWarnings(object@localData$precursorMz <-
                         as.numeric(res[, "precursor_mz_text"]))

    object@.tables <- list(
        msms_spectrum = colnames(
            dbGetQuery(dbcon, "select * from msms_spectrum limit 0")),
        ms_compound = colnames(
            dbGetQuery(dbcon, "select * from ms_compound limit 0")),
        synonym = colnames(
            dbGetQuery(dbcon, "select * from synonym")))
    object@spectraVariables <- c(.map_sql_to_spectraVariables(
        unique(unlist(object@.tables))),
        "precursor_mz_text", "compound_name")
    validObject(object)
    object
})

#' @exportMethod acquisitionNum
#'
#' @importMethodsFrom ProtGenerics acquisitionNum
#'
#' @rdname MsBackendMassbankSql
setMethod("acquisitionNum", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "acquisitionNum")[, 1]
    else integer()
})

#' @importMethodsFrom Spectra peaksData
#'
#' @exportMethod peaksData
#'
#' @rdname MsBackendMassbankSql
setMethod("peaksData", "MsBackendMassbankSql", function(object) {
    pks <- .fetch_peaks_sql(object)
    f <- factor(pks$spectrum_id)
    pks <- unname(split.data.frame(pks, f)[object@spectraIds])
    lapply(pks, function(z) {
        if (nrow(z))
            as.matrix(z[, 2:3], rownames.force = FALSE)
        else matrix(ncol = 2, nrow = 0,
                    dimnames = list(character(), c("mz", "intensity")))
    })
})

#' @exportMethod centroided
#'
#' @aliases centroided<-,MsBackendMassbankSql-method
#'
#' @importMethodsFrom ProtGenerics centroided
#'
#' @rdname MsBackendMassbankSql
setMethod("centroided", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "centroided")[, 1]
    else logical()
})

#' @exportMethod centroided<-
#'
#' @importMethodsFrom ProtGenerics centroided<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("centroided", "MsBackendMassbankSql", function(object, value) {
    if (!is.logical(value))
        stop("'value' has to be a logical")
    object$centroided <- value
    validObject(object)
    object
})

#' @exportMethod collisionEnergy
#'
#' @importMethodsFrom ProtGenerics collisionEnergy
#'
#' @rdname MsBackendMassbankSql
setMethod("collisionEnergy", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "collisionEnergy")[, 1]
    else numeric()
})

#' @exportMethod collisionEnergy<-
#'
#' @importMethodsFrom ProtGenerics collisionEnergy<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("collisionEnergy", "MsBackendMassbankSql",
                 function(object, value) {
                     if (!is.numeric(value))
                         stop("'value' has to be a numeric value")
                     object$collisionEnergy <- value
                     validObject(object)
                     object
})

#' @exportMethod dataOrigin
#'
#' @importMethodsFrom ProtGenerics dataOrigin
#'
#' @rdname MsBackendMassbankSql
setMethod("dataOrigin", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "dataOrigin")[, 1]
    else character()
})

#' @exportMethod dataOrigin<-
#'
#' @importMethodsFrom ProtGenerics dataOrigin<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("dataOrigin", "MsBackendMassbankSql", function(object, value) {
    if (!is.character(value))
        stop("'value' has to be a character")
    object$dataOrigin <- value
    validObject(object)
    object
})

#' @exportMethod dataStorage
#'
#' @importMethodsFrom ProtGenerics dataStorage
#'
#' @rdname MsBackendMassbankSql
setMethod("dataStorage", "MsBackendMassbankSql", function(object) {
    rep("<MassBank>", length(object))
})

#' @exportMethod intensity
#'
#' @importMethodsFrom ProtGenerics intensity
#'
#' @rdname MsBackendMassbankSql
setMethod("intensity", "MsBackendMassbankSql", function(object) {
    if (length(object)) {
        .spectra_data_massbank_sql(object, "intensity")[, 1]
    } else NumericList()
})

#' @exportMethod intensity<-
#'
#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("intensity", "MsBackendMassbankSql", function(object, value) {
    stop("Can not replace original intensity values in MassBank.")
})

#' @exportMethod ionCount
#'
#' @importMethodsFrom ProtGenerics ionCount
#'
#' @importFrom MsCoreUtils vapply1d
#'
#' @rdname MsBackendMassbankSql
setMethod("ionCount", "MsBackendMassbankSql", function(object) {
    vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @exportMethod isEmpty
#'
#' @rdname MsBackendMassbankSql
#'
#' @importMethodsFrom S4Vectors isEmpty
setMethod("isEmpty", "MsBackendMassbankSql", function(x) {
    lengths(intensity(x)) == 0
})

#' @exportMethod isolationWindowLowerMz
#'
#' @importMethodsFrom ProtGenerics isolationWindowLowerMz
#'
#' @rdname MsBackendMassbankSql
setMethod("isolationWindowLowerMz", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "isolationWindowLowerMz")[, 1]
    else numeric()
})

#' @exportMethod isolationWindowLowerMz<-
#'
#' @importMethodsFrom ProtGenerics isolationWindowLowerMz<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("isolationWindowLowerMz", "MsBackendMassbankSql",
                 function(object, value) {
                     if (!is.numeric(value))
                         stop("'value' has to be numeric")
                     object$isolationWindowLowerMz <- value
                     validObject(object)
                     object
                 })

#' @exportMethod isolationWindowTargetMz
#'
#' @importMethodsFrom ProtGenerics isolationWindowTargetMz
#'
#' @rdname MsBackendMassbankSql
setMethod("isolationWindowTargetMz", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "isolationWindowTargetMz")[, 1]
    else numeric()
})

#' @exportMethod isolationWindowTargetMz<-
#'
#' @importMethodsFrom ProtGenerics isolationWindowTargetMz<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("isolationWindowTargetMz", "MsBackendMassbankSql",
                 function(object, value) {
                     if (!is.numeric(value))
                         stop("'value' has to be numeric")
                     object$isolationWindowTargetMz <- value
                     validObject(object)
                     object
                 })

#' @exportMethod isolationWindowUpperMz
#'
#' @importMethodsFrom ProtGenerics isolationWindowUpperMz
#'
#' @rdname MsBackendMassbankSql
setMethod("isolationWindowUpperMz", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "isolationWindowUpperMz")[, 1]
    else numeric()
})

#' @exportMethod isolationWindowUpperMz<-
#'
#' @importMethodsFrom ProtGenerics isolationWindowUpperMz<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("isolationWindowUpperMz", "MsBackendMassbankSql",
                 function(object, value) {
                     if (!is.numeric(value))
                         stop("'value' has to be numeric")
                     object$isolationWindowUpperMz <- value
                     validObject(object)
                     object
                 })

#' @exportMethod length
#'
#' @rdname MsBackendMassbankSql
setMethod("length", "MsBackendMassbankSql", function(x) {
    length(x@spectraIds)
})

#' @exportMethod msLevel
#'
#' @importMethodsFrom ProtGenerics msLevel
#'
#' @rdname MsBackendMassbankSql
setMethod("msLevel", "MsBackendMassbankSql", function(object) {
    if (length(object)) {
        .spectra_data_massbank_sql(object, "msLevel")[, 1]
    } else integer()
})

#' @exportMethod mz
#'
#' @importMethodsFrom ProtGenerics mz
#'
#' @rdname MsBackendMassbankSql
setMethod("mz", "MsBackendMassbankSql", function(object) {
    if (length(object)) {
        .spectra_data_massbank_sql(object, "mz")[, 1]
    } else NumericList()
})

#' @exportMethod mz<-
#'
#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("mz", "MsBackendMassbankSql", function(object, value) {
    stop("Can not replace original data in MassBank.")
})

#' @rdname MsBackendMassbankSql
setMethod("lengths", "MsBackendMassbankSql", function(x, use.names = FALSE) {
    lengths(mz(x))
})

#' @exportMethod polarity
#'
#' @importMethodsFrom ProtGenerics polarity
#'
#' @rdname MsBackendMassbankSql
setMethod("polarity", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "polarity")[, 1]
    else integer()
})

#' @exportMethod polarity<-
#'
#' @importMethodsFrom ProtGenerics polarity<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("polarity", "MsBackendMassbankSql", function(object, value) {
    if (!is.numeric(value))
        stop("'value' has to be numeric")
    object$polarity <- as.integer(value)
    validObject(object)
    object
})

#' @exportMethod precursorCharge
#'
#' @importMethodsFrom ProtGenerics precursorCharge
#'
#' @rdname MsBackendMassbankSql
setMethod("precursorCharge", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "precursorCharge")[, 1]
    else integer()
})

#' @exportMethod precursorIntensity
#'
#' @importMethodsFrom ProtGenerics precursorIntensity
#'
#' @rdname MsBackendMassbankSql
setMethod("precursorIntensity", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "precursorIntensity")[, 1]
    else numeric()
})

#' @exportMethod precursorMz
#'
#' @importMethodsFrom ProtGenerics precursorMz
#'
#' @rdname MsBackendMassbankSql
setMethod("precursorMz", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "precursorMz")[, 1]
    else numeric()
})

#' @exportMethod reset
#'
#' @importMethodsFrom Spectra reset
#'
#' @rdname MsBackendMassbankSql
setMethod("reset", "MsBackendMassbankSql", function(object) {
    message("Restoring original data ...", appendLF = FALSE)
    object@localData <- DataFrame()
    if (is(object@dbcon, "DBIConnection"))
        object <- backendInitialize(object, object@dbcon)
    message("DONE")
    object
})

#' @exportMethod rtime
#'
#' @importMethodsFrom ProtGenerics rtime
#'
#' @rdname MsBackendMassbankSql
setMethod("rtime", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "rtime")[, 1]
    else numeric()
})

#' @exportMethod rtime<-
#'
#' @importMethodsFrom ProtGenerics rtime<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("rtime", "MsBackendMassbankSql", function(object, value) {
    if (!is.numeric(value))
        stop("'value' has to be numeric")
    object$rtime <- value
    validObject(object)
    object
})

#' @exportMethod scanIndex
#'
#' @importMethodsFrom ProtGenerics scanIndex
#'
#' @rdname MsBackendMassbankSql
setMethod("scanIndex", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "scanIndex")[, 1]
    else numeric()
})

#' @importMethodsFrom Spectra selectSpectraVariables
#'
#' @exportMethod selectSpectraVariables
#'
#' @rdname MsBackendMassbankSql
setMethod(
    "selectSpectraVariables", "MsBackendMassbankSql",
    function(object, spectraVariables = spectraVariables(object)) {
        if (any(!spectraVariables %in% spectraVariables(object)))
            stop("spectra variable(s) ",
                 paste(spectraVariables[!spectraVariables %in%
                                        spectraVariables(object)],
                       collapse = ", "), " not available")
        object@spectraVariables <- intersect(object@spectraVariables,
                                             spectraVariables)
        object@coreSpectraVariables <- intersect(object@coreSpectraVariables,
                                                 spectraVariables)
        object@localData <- object@localData[, colnames(object@localData) %in%
                                               spectraVariables, drop = FALSE]
        validObject(object)
        object
    })

#' @exportMethod smoothed
#'
#' @importMethodsFrom ProtGenerics smoothed
#'
#' @rdname MsBackendMassbankSql
setMethod("smoothed", "MsBackendMassbankSql", function(object) {
    if (length(object))
        .spectra_data_massbank_sql(object, "smoothed")[, 1]
    else logical()
})

#' @exportMethod smoothed<-
#'
#' @aliases smoothed<-,MsBackendMassbankSql-method
#'
#' @importMethodsFrom ProtGenerics smoothed<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("smoothed", "MsBackendMassbankSql", function(object, value) {
    if (!is.logical(value))
        stop("'value' has to be logical")
    object$smoothed <- value
    validObject(object)
    object
})

#' @exportMethod spectraData
#'
#' @rdname MsBackendMassbankSql
setMethod(
    "spectraData", "MsBackendMassbankSql",
    function(object, columns = spectraVariables(object)) {
        .spectra_data_massbank_sql(object, columns = columns)
    })

#' @exportMethod spectraData<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("spectraData", "MsBackendMassbankSql",function(object, value) {
    stop(class(object)[1], " does not support replacing the full spectra data.")
})

#' @exportMethod spectraNames
#'
#' @importMethodsFrom ProtGenerics spectraNames
#'
#' @rdname MsBackendMassbankSql
setMethod("spectraNames", "MsBackendMassbankSql", function(object) {
    object@spectraIds
})

#' @exportMethod spectraNames<-
#'
#' @importMethodsFrom ProtGenerics spectraNames<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("spectraNames", "MsBackendMassbankSql",
                 function(object, value) {
                     stop(class(object)[1],
                          " does not support replacing spectra names (IDs).")
})

#' @exportMethod spectraVariables
#'
#' @importMethodsFrom ProtGenerics spectraVariables
#'
#' @rdname MsBackendMassbankSql
setMethod("spectraVariables", "MsBackendMassbankSql", function(object) {
    unique(c(object@coreSpectraVariables, colnames(object@localData),
             object@spectraVariables))
})

#' @exportMethod tic
#'
#' @importMethodsFrom ProtGenerics tic
#'
#' @rdname MsBackendMassbankSql
setMethod("tic", "MsBackendMassbankSql", function(object, initial = TRUE) {
    if (initial) {
        if (any(colnames(object@localData) == "totIonCurrent"))
            object@localData[, "totIonCurrent"]
        else rep(NA_real_, times = length(object))
    } else vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @exportMethod [
#'
#' @importFrom MsCoreUtils i2index
#'
#' @importFrom methods slot<-
#'
#' @importFrom S4Vectors extractROWS
#'
#' @rdname MsBackendMassbankSql
setMethod("[", "MsBackendMassbankSql", function(x, i, j, ..., drop = FALSE) {
    if (missing(i))
        return(x)
    i <- i2index(i, length(x), x@spectraIds)
    slot(x, "spectraIds", check = FALSE) <- x@spectraIds[i]
    if (length(x@localData))
        slot(x, "localData", check = FALSE) <- extractROWS(x@localData, i)
    x
})

#' @exportMethod $
#'
#' @rdname MsBackendMassbankSql
setMethod("$", "MsBackendMassbankSql", function(x, name) {
    if (!any(spectraVariables(x) == name))
        stop("Spectra variable '", name, "' not available.")
    .spectra_data_massbank_sql(x, name)[, 1]
})

#' @exportMethod $<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("$", "MsBackendMassbankSql", function(x, name, value) {
    if (name %in% c("mz", "intensity"))
        stop("Replacing m/z and intensity values is not supported.")
    if (length(value) == 1)
        value <- rep(value, length(x))
    if (length(value) != length(x))
        stop("value has to be either of length 1 or length equal to the ",
             "number of spectra")
    if (length(x@localData)) {
        cn <- colnames(x@localData) == name
        if (any(cn))
            x@localData[, cn] <- value
        else {
            cn <- colnames(x@localData)
            x@localData <- cbind(x@localData, value)
            colnames(x@localData) <- c(cn, name)
        }
    } else {
        x@localData <- DataFrame(value)
        colnames(x@localData) <- name
    }
    validObject(x)
    x
})
