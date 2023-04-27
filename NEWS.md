# MsBackendSql 1.1

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
