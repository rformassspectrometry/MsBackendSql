# `Spectra` MS backend storing data in a SQL database

The `MsBackendSql` is an implementation for the
[`Spectra::MsBackend()`](https://rdrr.io/pkg/Spectra/man/MsBackend.html)
class for
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
objects which stores and retrieves MS data from a SQL database. New
databases can be created from raw MS data files using
`createMsBackendSqlDatabase()`.

## Usage

``` r
MsBackendSql()

createMsBackendSqlDatabase(
  dbcon,
  x = character(),
  backend = MsBackendMzR(),
  chunksize = 10L,
  blob = TRUE,
  peaksStorageMode = c("blob2", "long", "blob"),
  partitionBy = c("none", "spectrum", "chunk"),
  partitionNumber = 10L
)

# S4 method for class 'MsBackendSql'
show(object)

# S4 method for class 'MsBackendSql'
backendInitialize(object, dbcon, data, ...)

# S4 method for class 'MsBackendSql'
dataStorage(object)

# S4 method for class 'MsBackendSql'
x[i, j, ..., drop = FALSE]

# S4 method for class 'MsBackendSql,ANY'
extractByIndex(object, i)

# S4 method for class 'MsBackendSql'
peaksData(object, columns = c("mz", "intensity"))

# S4 method for class 'MsBackendSql'
peaksVariables(object)

# S4 method for class 'MsBackendSql'
intensity(object)

# S4 method for class 'MsBackendSql'
intensity(object) <- value

# S4 method for class 'MsBackendSql'
mz(object)

# S4 method for class 'MsBackendSql'
mz(object) <- value

# S4 method for class 'MsBackendSql'
x$name <- value

# S4 method for class 'MsBackendSql'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MsBackendSql'
reset(object)

# S4 method for class 'MsBackendSql'
spectraNames(object)

# S4 method for class 'MsBackendSql'
spectraNames(object) <- value

# S4 method for class 'MsBackendSql'
filterMsLevel(object, msLevel = uniqueMsLevels(object))

# S4 method for class 'MsBackendSql'
filterRt(object, rt = numeric(), msLevel. = integer())

# S4 method for class 'MsBackendSql'
filterDataOrigin(object, dataOrigin = character())

# S4 method for class 'MsBackendSql'
filterPrecursorMzRange(object, mz = numeric())

# S4 method for class 'MsBackendSql'
filterPrecursorMzValues(object, mz = numeric(), ppm = 20, tolerance = 0)

# S4 method for class 'MsBackendSql'
uniqueMsLevels(object, ...)

# S4 method for class 'MsBackendSql'
backendMerge(object, ...)

# S4 method for class 'MsBackendSql'
precScanNum(object)

# S4 method for class 'MsBackendSql'
centroided(object)

# S4 method for class 'MsBackendSql'
smoothed(object)

# S4 method for class 'MsBackendSql'
tic(object, initial = TRUE)

# S4 method for class 'MsBackendSql'
supportsSetBackend(object, ...)

# S4 method for class 'MsBackendSql'
backendBpparam(object, BPPARAM = bpparam())

# S4 method for class 'MsBackendSql'
dbconn(x)

# S4 method for class 'MsBackendSql'
longForm(object, columns = spectraVariables(object))
```

## Arguments

- dbcon:

  Connection to a database.

- x:

  For `createMsBackendSqlDatabase()`: `character` with the names of the
  raw data files from which the data should be imported. For other
  methods an `MsqlBackend` instance.

- backend:

  For `createMsBackendSqlDatabase()`: MS backend that can be used to
  import MS data from the raw files specified with parameter `x`.

- chunksize:

  For `createMsBackendSqlDatabase()`: `integer(1)` defining the number
  of input that should be processed per iteration. With `chunksize = 1`
  each file specified with `x` will be imported and its data inserted to
  the database. With `chunksize = 5` data from 5 files will be imported
  (in parallel) and inserted to the database. Thus, higher values might
  result in faster database creation, but require also more memory.

- blob:

  For `createMsBackendSqlDatabase()`,
  [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  `logical(1)` whether individual m/z and intensity values should be
  stored separately (`blob = FALSE`) or if the peaks data should be
  stored as data type *BLOB* into the database (`blob = TRUE`, the
  default). See also parameter `peaksStorageMode` for different data
  storage options.

- peaksStorageMode:

  `character(1)` defining how peaks variables are stored in the
  database. The default `peaksStorageMode = "blob2"` stores the full
  peaks matrix of each spectrum as data type *BLOB* as one entry into a
  single database table column. `peaksStorageMode = "blob"` stores in
  contrast the m/z and intensity vectors as separate BLOB types into two
  database tables. `peaksStorageMode = "long"` allows to store the data
  in *long mode*, i.e. each intensity and m/z value is stored
  individually in the database.

- partitionBy:

  For `createMsBackendSqlDatabase()`: `character(1)` defining if and how
  the peak data table should be partitioned. `"none"` (default): no
  partitioning, `"spectrum"`: peaks are assigned to the partition based
  on the spectrum ID (number), i.e. spectra are evenly (consecutively)
  assigned across partitions. For `partitionNumber = 3`, the first
  spectrum is assigned to the first partition, the second to the second,
  the third to the third and the fourth spectrum again to the first
  partition. `"chunk"`: spectra processed as part of the same *chunk*
  are placed into the same partition. All spectra from the next
  processed chunk are assigned to the next partition. Note that this is
  only available for MySQL/MariaDB databases, i.e., if `con` is a
  `MariaDBConnection`. See details for more information.

- partitionNumber:

  For `createMsBackendSqlDatabase()`: `integer(1)` defining the number
  of partitions the database table will be partitioned into (only
  supported for MySQL/MariaDB databases).

- object:

  A `MsBackendSql` instance.

- data:

  For
  [`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  optional `DataFrame` with the full spectra data that should be
  inserted into a (new) `MsBackendSql` database. If provided, it is
  assumed that `dbcon` is a (writeable) connection to an empty database
  into which `data` should be inserted. `data` could be the output of
  `spectraData` from another backend.

- ...:

  For `[`: ignored. For
  [`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html),
  if parameter `data` is used: additional parameters to be passed to the
  function creating the database such as `blob` or `peaksStorageMode`.
  For
  [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  any parameters supported by
  [`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  or `createMsBackendSqlDatabase()`.

- i:

  For `[`: `integer` or `logical` to subset the object.

- j:

  For `[`: ignored.

- drop:

  For `[`: `logical(1)`, ignored.

- columns:

  For
  [`spectraData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  [`character()`](https://rdrr.io/r/base/character.html) optionally
  defining a subset of spectra variables that should be returned.
  Defaults to `columns = spectraVariables(object)` hence all variables
  are returned. For `peaksData` accessor: optional `character` with
  requested columns in the individual `matrix` of the returned `list`.
  Defaults to `columns = c("mz", "intensity")` but all columns listed by
  `peaksVariables` would be supported.

- value:

  For all setter methods: replacement value.

- name:

  For `<-`: `character(1)` with the name of the spectra variable to
  replace.

- msLevel:

  For
  [`filterMsLevel()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `integer` specifying the MS levels to filter the data.

- rt:

  For
  [`filterRt()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `numeric(2)` with the lower and upper retention time. Spectra with a
  retention time `>= rt[1]` and `<= rt[2]` are returned.

- msLevel.:

  For
  `filterRt(): `integer`with the MS level(s) on which the retention time filter should be applied (all spectra from other MS levels are considered for the filter and are returned *as is*). If not specified, the retention time filter is applied to all MS levels in`object\`.

- dataOrigin:

  For
  [`filterDataOrigin()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `character` with *data origin* values to which the data should be
  subsetted.

- mz:

  For
  [`filterPrecursorMzRange()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `numeric(2)` with the desired lower and upper limit of the precursor
  m/z range. For
  [`filterPrecursorMzValues()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `numeric` with the m/z value(s) to filter the object.

- ppm:

  For
  [`filterPrecursorMzValues()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `numeric` with the m/z-relative maximal acceptable difference for a
  m/z value to be considered matching. Can be of length 1 or equal to
  `length(mz)`.

- tolerance:

  For
  [`filterPrecursorMzValues()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `numeric` with the absolute difference for m/z values to be considered
  matching. Can be of length 1 or equal to `length(mz)`.

- initial:

  For [`tic()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `logical(1)` whether the original total ion count should be returned
  (`initial = TRUE`, the default) or whether it should be calculated on
  the spectras' intensities (`initial = FALSE`).

- BPPARAM:

  for
  [`backendBpparam()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  `BiocParallel` parallel processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

## Value

See documentation of respective function.

## Details

The `MsBackendSql` class is principally a *read-only* backend but by
extending the
[`Spectra::MsBackendCached()`](https://rdrr.io/pkg/Spectra/man/MsBackendCached.html)
backend from the `Spectra` package it allows changing and adding
(**temporarily**) spectra variables **without** changing the original
data in the SQL database.

## Note

The `MsBackendSql` backend keeps an (open) connection to the SQL
database with the data and hence does not support saving/loading of a
backend to disk (e.g. using `save` or `saveRDS`). Also, for the same
reason, the `MsBackendSql` does not support parallel processing. The
[`backendBpparam()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method for `MsBackendSql` will thus always return a
[`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
object.

The
[`MsBackendOfflineSql()`](https://rformassspectrometry.github.io/MsBackendSql/reference/MsBackendOfflineSql.md)
could be used as an alternative as it supports saving/loading the data
to/from disk and supports also parallel processing.

## Creation of backend objects

New backend objects can be created with the `MsBackendSql()` function.
SQL databases can be created and filled with MS data from raw data files
using the `createMsBackendSqlDatabase()` function or using
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
and providing all data with parameter `data`. In addition it is possible
to create a database from a `Spectra` object changing its backend to a
`MsBackendSql` or `MsBackendOfflineSql` using the
[`Spectra::setBackend()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
function. Existing SQL databases (created previously with
`createMsBackendSqlDatabase()` or
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
with the `data` parameter) can be loaded using the *conventional* way to
create/initialize `MsBackend` classes, i.e. using
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html).

- `createMsBackendSqlDatabase()`: create a database and fill it with MS
  data. Parameter `dbcon` is expected to be a database connection,
  parameter `x` a `character` vector with the file names from which to
  import the data. Parameter `backend` is used for the actual data
  import and defaults to `backend = MsBackendMzR()` hence allowing to
  import data from mzML, mzXML or netCDF files. Parameter `chunksize`
  allows to define the number of files (`x`) from which the data should
  be imported in one iteration. With the default `chunksize = 10L` data
  is imported from 10 files in `x` at the same time (if `backend`
  supports it even in parallel) and this data is then inserted into the
  database. Larger chunk sizes will require more memory and also larger
  disk space (as data import is performed through temporary files) but
  might eventually be faster. Parameter `blob` allows to define whether
  m/z and intensity values from a spectrum should be stored as a *BLOB*
  SQL data type in the database (`blob = TRUE`, the default) or if
  individual m/z and intensity values for each peak should be stored
  separately (`blob = FALSE`). The latter case results in a much larger
  database and slower performance of the `peaksData` function, but would
  allow to define custom (manual) SQL queries on individual peak values.
  For `blob = TRUE`, the peaks data can be stored in two different ways
  which can be selected with the additional parameter
  `peaksStorageMode`. The default `peaksStorageMode = "blob2"` stores
  the full peaks matrix of a spectrum as a single entry to the database
  (into a single database table column) while
  `peaksStorageMode = "blob"` (which was the default until version
  1.7.2) stores the m/z and intensity vectors as BLOB data type into two
  separate database table columns. Performance for
  [`peaksData()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html)
  is thus about twice as fast for `peaksStorageMode = "blob2"`. While
  data can be stored in any SQL database, at present it is suggested to
  use MySQL/MariaDB databases. For `dbcon` being a connection to a
  MySQL/MariaDB database, the tables will use the *ARIA* engine
  providing faster data access and will use *table partitioning*: tables
  are splitted into multiple partitions which can improve data insertion
  and index generation. Partitioning can be defined with the parameters
  `partitionBy` and `partitionNumber`. By default `partitionBy = "none"`
  no partitioning is performed. For `blob = TRUE` partitioning is
  usually not required. Only for `blob = FALSE ` and very large datasets
  it is suggested to enable table partitioning by selecting either
  `partitionBy = "spectrum"` or `partitionBy = "chunk"`. The first
  option assignes consecutive spectra to different partitions while the
  latter puts spectra from files part of the same *chunk* into the same
  partition. Both options have about the same performance but
  `partitionBy = "spectrum"` requires less disk space. Note that, while
  inserting the data takes a considerable amount of time, also the
  subsequent creation of database indices can take very long (even
  longer than data insertion for `blob = FALSE`).

- [`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  get access and initialize a `MsBackendSql` object. Parameter `object`
  is supposed to be a `MsBackendSql` instance, created e.g. with
  `MsBackendSql()`. Parameter `dbcon` is expected to be a connection to
  an existing *MsBackendSql* SQL database (created e.g. with
  `createMsBackendSqlDatabase()`). To initialize a
  [`MsBackendOfflineSql()`](https://rformassspectrometry.github.io/MsBackendSql/reference/MsBackendOfflineSql.md)
  all information required for a
  [`DBI::dbConnect()`](https://dbi.r-dbi.org/reference/dbConnect.html)
  call to connect to a database need to be provided.
  [`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  can alternatively also be used to create a **new** `MsBackendSql`
  database using the optional `data` parameter. In this case, `dbcon` is
  expected to be a writeable connection to an empty database and `data`
  a `DataFrame` with the **full** spectra and peaks data to be inserted
  into this database. The format of `data` should match the format of
  the `DataFrame` returned by the
  [`spectraData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  function and requires columns `"mz"` and `"intensity"` with the m/z
  and intensity values of each spectrum. The
  [`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  call will then create all necessary tables in the database, will fill
  these tables with the provided data and will return an `MsBackendSql`
  for this database. The `MsBackendSql` and `MsBackendOfflineSql`
  objects support the
  [`Spectra::setBackend()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  method from `Spectra` to change from (any) backend to a
  `MsBackendSql`. Any parameters to the
  [`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  function can be passed to
  [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html).
  Note however that chunk-wise (or parallel) processing needs to be
  disabled in this case by passing eventually `f = factor()` to the
  [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  call.

- [`supportsSetBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  whether `MsBackendSql` supports the
  [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  method to change the `MsBackend` of a `Spectra` object to a
  `MsBackendSql`. Returns `TRUE`, thus, changing the backend to a
  `MsBackendSql` is supported **if** a writeable database connection is
  provided in addition with parameter `dbcon` (i.e.
  `setBackend(sps, MsBackendSql(), dbcon = con)` with `con` being a
  connection to an **empty** database would store the full spectra data
  from the `Spectra` object `sps` into the specified database and would
  return a `Spectra` object that uses a `MsBackendSql`).

- [`backendBpparam()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  whether a `MsBackendSql` supports parallel processing. Takes a
  `MsBackendSql` and a parallel processing setup (see
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for details) as input and always returns a
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  since `MsBackendSql` does **not** support parallel processing.

- `dbconn()`: returns the connection to the database.

## Subsetting, merging and filtering data

`MsBackendSql` objects can be subsetted using the `[` or
[`extractByIndex()`](https://rdrr.io/pkg/ProtGenerics/man/extractByIndex.html)
functions. Internally, this will simply subset the `integer` vector of
the primary keys and eventually cached data. The original data in the
database **is not** affected by any subsetting operation. Any subsetting
operation can be *undone* by resetting the object with the
[`reset()`](https://rdrr.io/pkg/Spectra/man/addProcessing.html)
function. Subsetting in arbitrary order as well as index replication is
supported.

Multiple `MsBackendSql` objects can also be merged (combined) with the
`backendMerge()` function. Note that this requires that all
`MsBackendSql` objects are connected to the **same** database. This
function is thus mostly used for combining `MsBackendSql` objects that
were previously splitted using e.g.
[`split()`](https://rdrr.io/r/base/split.html).

In addition, `MsBackendSql` supports all other filtering methods
available through
[`Spectra::MsBackendCached()`](https://rdrr.io/pkg/Spectra/man/MsBackendCached.html).
Implementation of filter functions optimized for `MsBackendSql` objects
are:

- [`filterDataOrigin()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  filter the object retaining spectra with `dataOrigin` spectra variable
  values matching the provided ones with parameter `dataOrigin`. The
  function returns the results in the order of the values provided with
  parameter `dataOrigin`.

- [`filterMsLevel()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  filter the object based on the MS levels specified with parameter
  `msLevel`. The function does the filtering using SQL queries. If
  `"msLevel"` is a *local* variable stored within the object (and hence
  in memory) the default implementation in `MsBackendCached` is used
  instead.

- [`filterPrecursorMzRange()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  filters the data keeping only spectra with a `precursorMz` within the
  m/z value range provided with parameter `mz` (i.e. all spectra with a
  precursor m/z `>= mz[1L]` and `<= mz[2L]`).

- filterPrecursorMzValues()`: filters the data keeping only spectra with precursor m/z values matching the value(s) provided with parameter `mz`. Parameters `ppm`and`tolerance`allow to specify acceptable differences between compared values. Lengths of`ppm`and`tolerance`can be either`1`or equal to`length(mz)\`
  to use different values for ppm and tolerance for each provided m/z
  value.

- [`filterRt()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  filter the object keeping only spectra with retention times within the
  specified retention time range (parameter `rt`). Optional parameter
  `msLevel.` allows to restrict the retention time filter only on the
  provided MS level(s) returning all spectra from other MS levels.

## Accessing and *modifying* data

The functions listed here are specifically implemented for
`MsBackendSql`. In addition, `MsBackendSql` inherits and supports all
data accessor, filtering functions and data manipulation functions from
[`Spectra::MsBackendCached()`](https://rdrr.io/pkg/Spectra/man/MsBackendCached.html).

- `$`, `$<-`: access or set (add) spectra variables in `object`. Spectra
  variables added or modified using the `$<-` are *cached* locally
  within the object (data in the database is never changed). To restore
  an object (i.e. drop all cached values) the `reset` function can be
  used.

- [`dataStorage()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  returns a `character` vector same length as there are spectra in
  `object` with the name of the database containing the data.

- `intensity<-`: not supported.

- `longForm()`: extract the MS data in *long form* as a `data.frame`.
  Parameter `columns` allows to specify the columns (spectra and/or
  peaks variables) that should be included in the result. If MS peaks
  data are stored in long form in the database (i.e.,
  `peaksStorageMode = "long"` is used), the data is extracted using a
  dedicated SQL query. Otherwise the default implementation of the
  `longForm()` method from the *Spectra* package that is based on the
  [`spectraData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  function is used instead. Note that the performance of the SQL-based
  `longForm()` function is not necessarily higher than the default
  implementation, mostly because data extraction from the database
  layouts that store the MS peaks data as *BLOB* datatype (i.e.,
  `peaksStorageMode = "blob"` or `peaksStorageMode = "blob2"`) is
  faster.

- `mz<-`: not supported.

- [`peaksData()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html):
  returns a `list` with the spectras' peak data. The length of the list
  is equal to the number of spectra in `object`. Each element of the
  list is a `matrix` with columns according to parameter `columns`. For
  an empty spectrum, a `matrix` with 0 rows is returned. Use
  `peaksVariables(object)` to list supported values for parameter
  `columns`.

- [`peaksVariables()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html):
  returns a `character` with the available peak variables, i.e. columns
  that could be queried with
  [`peaksData()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html).

- [`reset()`](https://rdrr.io/pkg/Spectra/man/addProcessing.html):
  *restores* an `MsBackendSql` by re-initializing it with the data from
  the database. Any subsetting or cached spectra variables will be lost.

- [`spectraData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  gets general spectrum metadata.
  [`spectraData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  returns a `DataFrame` with the same number of rows as there are
  spectra in `object`. Parameter `columns` allows to select specific
  spectra variables.

- [`spectraNames()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html),
  `spectraNames<-`: returns a `character` of length equal to the number
  of spectra in `object` with the primary keys of the spectra from the
  database (converted to `character`). Replacing spectra names with
  `spectraNames<-` is not supported.

- [`uniqueMsLevels()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  returns the unique MS levels of all spectra in `object`.

- [`tic()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  returns the originally reported total ion count (for `initial = TRUE`)
  or calculates the total ion count from the intensities of each
  spectrum (for `initial = FALSE`).

## Implementation notes

Internally, the `MsBackendSql` class contains only the primary keys for
all spectra stored in the SQL database. Keeping only these `integer` in
memory guarantees a minimal memory footpring of the object. Still,
depending of the number of spectra in the database, this `integer`
vector might become very large. Any data access will involve SQL calls
to retrieve the data from the database. By extending the
[`Spectra::MsBackendCached()`](https://rdrr.io/pkg/Spectra/man/MsBackendCached.html)
object from the `Spectra` package, the `MsBackendSql` supports to
(temporarily, i.e. for the duration of the R session) add or modify
spectra variables. These are however stored in a `data.frame` within the
object thus increasing the memory demand of the object.

## Author

Johannes Rainer

## Examples

``` r

####
## Create a new MsBackendSql database

## Define a file from which to import the data
data_file <- MsDataHub::X20171016_POOL_POS_3_105.134.mzML()
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache

## Create a database/connection to a database
library(RSQLite)
db_file <- tempfile()
dbc <- dbConnect(SQLite(), db_file)

## Import the data from the file into the database
createMsBackendSqlDatabase(dbc, data_file)
#> Importing data ... 
#> 
#> [==========================================================] 1/1 (100%) in  2s
#> 
#> Creating indices 
#> .
#> .
#> .
#> .
#>  Done
#> [1] TRUE
dbDisconnect(dbc)

## Initialize a MsBackendSql
dbc <- dbConnect(SQLite(), db_file)
be <- backendInitialize(MsBackendSql(), dbc)

be
#> MsBackendSql with 931 spectra
#>       msLevel precursorMz  polarity
#>     <integer>   <numeric> <integer>
#> 1           1          NA         1
#> 2           1          NA         1
#> 3           1          NA         1
#> 4           1          NA         1
#> 5           1          NA         1
#> ...       ...         ...       ...
#> 927         1          NA         1
#> 928         1          NA         1
#> 929         1          NA         1
#> 930         1          NA         1
#> 931         1          NA         1
#>  ... 35 more variables/columns.
#>  Use  'spectraVariables' to list all of them.
#> Database: /tmp/RtmpDEx8im/file16b77de1154e

## Original data source
head(be$dataOrigin)
#> [1] "/github/home/.cache/R/ExperimentHub/102835a6d416_7860"
#> [2] "/github/home/.cache/R/ExperimentHub/102835a6d416_7860"
#> [3] "/github/home/.cache/R/ExperimentHub/102835a6d416_7860"
#> [4] "/github/home/.cache/R/ExperimentHub/102835a6d416_7860"
#> [5] "/github/home/.cache/R/ExperimentHub/102835a6d416_7860"
#> [6] "/github/home/.cache/R/ExperimentHub/102835a6d416_7860"

## Data storage
head(dataStorage(be))
#> [1] "/tmp/RtmpDEx8im/file16b77de1154e" "/tmp/RtmpDEx8im/file16b77de1154e"
#> [3] "/tmp/RtmpDEx8im/file16b77de1154e" "/tmp/RtmpDEx8im/file16b77de1154e"
#> [5] "/tmp/RtmpDEx8im/file16b77de1154e" "/tmp/RtmpDEx8im/file16b77de1154e"

## Access all spectra data
spd <- spectraData(be)
spd
#> DataFrame with 931 rows and 38 columns
#>       msLevel     rtime acquisitionNum scanIndex                          mz
#>     <integer> <numeric>      <integer> <integer>               <NumericList>
#> 1           1     0.275              1         1 108.078,108.080,108.081,...
#> 2           1     0.554              2         2 105.040,105.041,105.043,...
#> 3           1     0.833              3         3 105.037,105.038,105.040,...
#> 4           1     1.112              4         4 105.037,105.038,105.040,...
#> 5           1     1.391              5         5 105.034,105.035,105.037,...
#> ...       ...       ...            ...       ...                         ...
#> 927         1   258.636            927       927 105.008,105.009,105.010,...
#> 928         1   258.915            928       928 105.009,105.010,105.012,...
#> 929         1   259.194            929       929 105.010,105.012,105.013,...
#> 930         1   259.473            930       930 105.012,105.013,105.015,...
#> 931         1   259.752            931       931 105.012,105.013,105.015,...
#>           intensity            dataStorage             dataOrigin centroided
#>       <NumericList>            <character>            <character>  <logical>
#> 1     0,412,  0,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#> 2     0,110,110,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#> 3     0,156,156,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#> 4     0,127,  0,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#> 5     0,123,123,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#> ...             ...                    ...                    ...        ...
#> 927     0,22,67,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#> 928     0,21, 0,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#> 929     0,26, 0,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#> 930     0,23, 0,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#> 931     0,25, 0,... /tmp/RtmpDEx8im/file.. /github/home/.cache/..      FALSE
#>      smoothed  polarity precScanNum precursorMz precursorIntensity
#>     <logical> <integer>   <integer>   <numeric>          <numeric>
#> 1          NA         1          NA          NA                 NA
#> 2          NA         1          NA          NA                 NA
#> 3          NA         1          NA          NA                 NA
#> 4          NA         1          NA          NA                 NA
#> 5          NA         1          NA          NA                 NA
#> ...       ...       ...         ...         ...                ...
#> 927        NA         1          NA          NA                 NA
#> 928        NA         1          NA          NA                 NA
#> 929        NA         1          NA          NA                 NA
#> 930        NA         1          NA          NA                 NA
#> 931        NA         1          NA          NA                 NA
#>     precursorCharge collisionEnergy isolationWindowLowerMz
#>           <integer>       <numeric>              <numeric>
#> 1                NA              NA                     NA
#> 2                NA              NA                     NA
#> 3                NA              NA                     NA
#> 4                NA              NA                     NA
#> 5                NA              NA                     NA
#> ...             ...             ...                    ...
#> 927              NA              NA                     NA
#> 928              NA              NA                     NA
#> 929              NA              NA                     NA
#> 930              NA              NA                     NA
#> 931              NA              NA                     NA
#>     isolationWindowTargetMz isolationWindowUpperMz peaksCount totIonCurrent
#>                   <numeric>              <numeric>  <integer>     <numeric>
#> 1                        NA                     NA        542        725013
#> 2                        NA                     NA       1777       1115422
#> 3                        NA                     NA       1360        979255
#> 4                        NA                     NA       1626        992255
#> 5                        NA                     NA       1694       1018768
#> ...                     ...                    ...        ...           ...
#> 927                      NA                     NA       3528        418159
#> 928                      NA                     NA       3662        429930
#> 929                      NA                     NA       3442        413427
#> 930                      NA                     NA       3564        447033
#> 931                      NA                     NA       3517        441025
#>     basePeakMZ basePeakIntensity electronBeamEnergy ionisationEnergy     lowMZ
#>      <numeric>         <numeric>          <numeric>        <numeric> <numeric>
#> 1      124.087            127807                 NA                0   108.078
#> 2      124.086            202258                 NA                0   105.040
#> 3      124.086            166954                 NA                0   105.037
#> 4      124.086            163449                 NA                0   105.037
#> 5      124.086            170723                 NA                0   105.034
#> ...        ...               ...                ...              ...       ...
#> 927    123.090             28544                 NA                0   105.008
#> 928    123.090             29164                 NA                0   105.009
#> 929    123.090             28724                 NA                0   105.010
#> 930    123.091             28743                 NA                0   105.012
#> 931    123.091             27503                 NA                0   105.012
#>        highMZ mergedScan mergedResultScanNum mergedResultStartScanNum
#>     <numeric>  <integer>           <integer>                <integer>
#> 1     133.988         NA                  NA                       NA
#> 2     133.982         NA                  NA                       NA
#> 3     133.987         NA                  NA                       NA
#> 4     133.987         NA                  NA                       NA
#> 5     133.992         NA                  NA                       NA
#> ...       ...        ...                 ...                      ...
#> 927   134.000         NA                  NA                       NA
#> 928   134.000         NA                  NA                       NA
#> 929   134.000         NA                  NA                       NA
#> 930   133.998         NA                  NA                       NA
#> 931   133.997         NA                  NA                       NA
#>     mergedResultEndScanNum injectionTime filterString             spectrumId
#>                  <integer>     <numeric>  <character>            <character>
#> 1                       NA             0           NA sample=1 period=1 cy..
#> 2                       NA             0           NA sample=1 period=1 cy..
#> 3                       NA             0           NA sample=1 period=1 cy..
#> 4                       NA             0           NA sample=1 period=1 cy..
#> 5                       NA             0           NA sample=1 period=1 cy..
#> ...                    ...           ...          ...                    ...
#> 927                     NA             0           NA sample=1 period=1 cy..
#> 928                     NA             0           NA sample=1 period=1 cy..
#> 929                     NA             0           NA sample=1 period=1 cy..
#> 930                     NA             0           NA sample=1 period=1 cy..
#> 931                     NA             0           NA sample=1 period=1 cy..
#>     ionMobilityDriftTime scanWindowLowerLimit scanWindowUpperLimit spectrum_id_
#>                <numeric>            <numeric>            <numeric>    <integer>
#> 1                     NA                   NA                   NA            1
#> 2                     NA                   NA                   NA            2
#> 3                     NA                   NA                   NA            3
#> 4                     NA                   NA                   NA            4
#> 5                     NA                   NA                   NA            5
#> ...                  ...                  ...                  ...          ...
#> 927                   NA                   NA                   NA          927
#> 928                   NA                   NA                   NA          928
#> 929                   NA                   NA                   NA          929
#> 930                   NA                   NA                   NA          930
#> 931                   NA                   NA                   NA          931

## Available variables
spectraVariables(be)
#>  [1] "msLevel"                  "rtime"                   
#>  [3] "acquisitionNum"           "scanIndex"               
#>  [5] "mz"                       "intensity"               
#>  [7] "dataStorage"              "dataOrigin"              
#>  [9] "centroided"               "smoothed"                
#> [11] "polarity"                 "precScanNum"             
#> [13] "precursorMz"              "precursorIntensity"      
#> [15] "precursorCharge"          "collisionEnergy"         
#> [17] "isolationWindowLowerMz"   "isolationWindowTargetMz" 
#> [19] "isolationWindowUpperMz"   "peaksCount"              
#> [21] "totIonCurrent"            "basePeakMZ"              
#> [23] "basePeakIntensity"        "electronBeamEnergy"      
#> [25] "ionisationEnergy"         "lowMZ"                   
#> [27] "highMZ"                   "mergedScan"              
#> [29] "mergedResultScanNum"      "mergedResultStartScanNum"
#> [31] "mergedResultEndScanNum"   "injectionTime"           
#> [33] "filterString"             "spectrumId"              
#> [35] "ionMobilityDriftTime"     "scanWindowLowerLimit"    
#> [37] "scanWindowUpperLimit"     "spectrum_id_"            

## Access mz values
mz(be)
#> NumericList of length 931
#> [[1]] 108.078350457464 108.079816722562 ... 133.985926725175 133.987559287009
#> [[2]] 105.039628321505 105.04107382602 ... 133.980645512489 133.982278040861
#> [[3]] 105.036737362206 105.038182846829 ... 133.985543127443 133.987175685653
#> [[4]] 105.036737362206 105.038182846829 ... 133.985543127443 133.987175685653
#> [[5]] 105.033846442691 105.035291907422 ... 133.990440831911 133.992073419959
#> [[6]] 105.033846442691 105.035291907422 ... 133.983910579179 133.985543127443
#> [[7]] 105.033846442691 105.035291907422 ... 133.982278040861 133.983910579179
#> [[8]] 105.036737362206 105.038182846829 ... 133.988808253809 133.990440831911
#> [[9]] 105.036737362206 105.038182846829 ... 133.985543127443 133.987175685653
#> [[10]] 105.039628321505 105.04107382602 ... 133.993706017952 133.995338625892
#> ...
#> <921 more elements>

## Subset the object to spectra in arbitrary order
be_sub <- be[c(5, 1, 1, 2, 4, 100)]
be_sub
#> MsBackendSql with 6 spectra
#>     msLevel precursorMz  polarity
#>   <integer>   <numeric> <integer>
#> 1         1          NA         1
#> 2         1          NA         1
#> 3         1          NA         1
#> 4         1          NA         1
#> 5         1          NA         1
#> 6         1          NA         1
#>  ... 35 more variables/columns.
#>  Use  'spectraVariables' to list all of them.
#> Database: /tmp/RtmpDEx8im/file16b77de1154e

## The internal spectrum IDs (primary keys from the database)
be_sub$spectrum_id_
#> [1]   5   1   1   2   4 100

## Add additional spectra variables
be_sub$new_variable <- "B"

## This variable is *cached* locally within the object (not inserted into
## the database)
be_sub$new_variable
#> [1] "B" "B" "B" "B" "B" "B"
```
