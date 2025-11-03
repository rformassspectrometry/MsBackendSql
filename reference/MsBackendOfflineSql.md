# SQL-based MS backend without active database connection

The `MsBackendOfflineSql` backend extends the
[`MsBackendSql()`](https://rformassspectrometry.github.io/MsBackendSql/reference/MsBackendSql.md)
backend directly and inherits thus all of its functions as well as
properties. The only difference between the two backend is that
`MsBackendSql` keeps an active connection to the SQL database inside the
object while the `MsBackendOfflineSql` backends reconnects to the SQL
database for each query. While the performance of the latter is slightly
lower (due to the need to connect/disconnect to the database for each
function call) it can also be used in a parallel processing environment.

## Usage

``` r
MsBackendOfflineSql()

# S4 method for class 'MsBackendOfflineSql'
backendInitialize(
  object,
  drv = NULL,
  dbname = character(),
  user = character(),
  password = character(),
  host = character(),
  port = NA_integer_,
  data,
  ...
)
```

## Arguments

- object:

  A `MsBackendOfflineSql` object.

- drv:

  A *DBI* database driver object (such as `SQLite()` from the `RSQLite`
  package or `MariaDB()` from the `RMariaDB` package). See
  [`DBI::dbConnect()`](https://dbi.r-dbi.org/reference/dbConnect.html)
  for more information.

- dbname:

  `character(1)` with the name of the database. Passed directly to
  [`DBI::dbConnect()`](https://dbi.r-dbi.org/reference/dbConnect.html).

- user:

  `character(1)` with the user name for the database. Passed directly to
  [`DBI::dbConnect()`](https://dbi.r-dbi.org/reference/dbConnect.html).

- password:

  `character(1)` with the password for the database. Note that this
  password is stored (unencrypted) within the object. Passed directly to
  [`DBI::dbConnect()`](https://dbi.r-dbi.org/reference/dbConnect.html).

- host:

  `character(1)` with the host running the database. Passed directly to
  [`DBI::dbConnect()`](https://dbi.r-dbi.org/reference/dbConnect.html).

- port:

  `integer(1)` with the port number (optional). Passed directly to
  [`DBI::dbConnect()`](https://dbi.r-dbi.org/reference/dbConnect.html).

- data:

  For
  [`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html):
  optional `DataFrame` with the full spectra data that should be
  inserted into a (new) `MsBackendSql` database. If provided, it is
  assumed that the provided database connection information if for a
  (writeable) empty database into which `data` should be inserted.
  `data` could be the output of `spectraData` from another backend.

- ...:

  ignored.

## Creation of backend objects

An empty instance of an `MsBackendOfflineSql` class can be created using
the `MsBackendOfflineSql()` function. An existing *MsBackendSql* SQL
database can be loaded with the
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
function. This function takes parameters `drv`, `dbname`, `user`,
`password`, `host` and `port`, all parameters that are passed to the
`dbConnect()` function to connect to the (**existing**) SQL database.

See
[`MsBackendSql()`](https://rformassspectrometry.github.io/MsBackendSql/reference/MsBackendSql.md)
for information on how to create a *MsBackend* SQL database.

## Author

Johannes Rainer
