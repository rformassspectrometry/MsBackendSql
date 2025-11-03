# MsBackendSql 1.10.1

- For storage mode of peaks data in long form (i.e., one row per peak), create a
  incremental (unique) *peak_id_* database column if the database requires that
  (e.g. for DuckDb).
- `longForm()` on a database requiring *peak_id_* (e.g. DuckDb) order results
  based on the *peak_id_* column.

# MsBackendSql 1.9

## Changes in 1.9.4

- Add additional unit tests to check results are similar to the reference
  implementation.

## Changes in 1.9.3

- Small performance improvement in `longForm()` implementation: avoid using
  `findMatches()` to order the result if not needed (i.e., if duplicated spectra
  are requested).

## Changes in 1.9.2

- Add a `longForm()` method for `MsBackendSql` to use an SQL query to extract
  the data in long from the database.

## Changes in 1.9.1

- Add support for *blob2* storage mode for *duckdb* databases.

# MsBackendSql 1.7

## Changes in 1.7.4

- Fix handling of parameter `msLevel.` in `filterRt()`.

## Changes in 1.7.3

- Add new peaks data storage mode *blob2* which stores the full peaks matrix
  as a single entity to the database table.
- Change default peaks data storage mode to *blob2*.
- Add parameter `peaksStorageMode` to database creation function allowing to
  select the new peaks data storage mode.
- Add `mz()` and `intensity()` methods.
- Small performance improvements for `peaksData()` using `fmatch()` from the
  *fastmatch* package and avoiding to re-order the results if not needed.
- Performance improvement for *blob* and *blob2* storage modes
  by using `xda = FALSE` in `serialize()`.
- Complete unit test coverage.

## Changes in 1.7.2

- Import `extractByIndex` from *ProtGenerics*.
- Add missing links in documentation.

## Changes in 1.7.1

- Complete unit tests to cover all code lines.

# MsBackendSql 1.5

## Changes in 1.5.1

- Implement the new `extractByIndex()` methods.

# MsBackendSql 1.3

## Changes in 1.3.5

- Improve input argument check and error message for `backendInitialize()` for
  `MsBackendOfflineSql`.
- Update documentation adding `()` to all function names.

## Changes in 1.3.4

- Ensure primary keys from the database are in the correct order for
  `backendInitialize()`.

## Changes in 1.3.3

- Import method generics from `ProtGenerics`.

## Changes in 1.3.2

- Add a dedicated `setBackend` method for `MsBackendSql` and
  `MsBackendOfflineSql` backends (issue
  [#17](https://github.com/rformassspectrometry/MsBackendSql/issues/17)).

## Changes in 1.3.1

- Add description on the use/advantages of different SQL database systems to the
  vignette.

# MsBackendSql 1.1

## Changes in 1.1.5

- Add `dbconn` methods for `MsBackendSql` and `MsBackendOfflineSql`.

## Changes in 1.1.4

- Improve performance of `createMsBackendSqlDatabase` by using also parallel
  processing for the `peaksData` call.

## Changes in 1.1.3

- Add support for `setBackend` to `MsBackendOfflineSql`.

## Changes in 1.1.2

- Mention in documentation that `MsBackendSql` can not be saved to disk.
- Expand vignette adding related documentation.

## Changes in 1.1.1

- Fix for `filterRt` avoiding to filter if range is infinite.


# MsBackendSql 0.99

## Changes in 0.99.7

- Decrease required R version to 4.2.

## Changes in 0.99.6

- Add `mzR` to *Suggests* to ensure package vignettes can be build properly.

## Changes in 0.99.5

- Add `MsBackendOfflineSql` backend that re-connects to the database for each
  query.

## Changes in 0.99.4

- Add `backendBpparam` method to ensure parallel processing is disabled for the
 `MsBackendSql` backend.

## Changes in 0.99.3

- Add `backendMerge` method.
- Add parameter `data` to `backendInitialize` to allow creating a new
  `MsBackendSql` database and store the values from `data` in it. This
  enables the use of `Spectra,setBackend` to convert any backend to a
  `MsBackendSql`.
- Implement `supportsSetBackend` to enable `setBackend,Spectra,MsBackendSql`.

## Changes in 0.99.2

- Evaluate validity of the `MsBackendSql` using the full unit test suite from
  the `Spectra` package.

## Changes in 0.99.1

- Address Kayla's package review comments.

## Changes in 0.99.0

- Prepare the package for submission to Bioconductor.


# MsBackendSql 0.98

## Changes in 0.98.1

- Rename `MsqlBackend` to `MsBackendSql`.

## Changes in 0.98.0

- Add vignette.


# MsBackendSql 0.1

## Changes in 0.1.0

- Add parameter `blob` to allow storing m/z and intensity values as *BLOB* data
  type in the database. `MsBackendSql` will use different functions to retrieve
  data from a database with this type of storage.

# MsBackendSql 0.0

## Changes in 0.0.5

- Add `filterMsLevel` method.
- Add `filterRt` method.
- Add `filterDataOrigin` method.
- Add `filterPrecursorMzRange` method.
- Add `filterPrecursorMzValues` method.

## Changes in 0.0.4

- Add optional parameter `tempfile`.

## Changes in 0.0.3

- Add `peaksVariables` methods and add support for parameter `columns` in
  `peaksData` (Requires `Spectra` version 1.15.17).

## Changes in 0.0.2

- First full implementation of `MsBackendSql`.
