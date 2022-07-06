# MsqlBackend 0.1

## Changes in 0.1.0

- Add parameter `blob` to allow storing m/z and intensity values as *BLOB* data
  type in the database. `MsqlBackend` will use different functions to retrieve
  data from a database with this type of storage.

# MsqlBackend 0.0

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

- First full implementation of `MsqlBackend`.
