# MsBackendSql 0.99

## Changes in 0.99.

- Add `backendMerge` method.

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
