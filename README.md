# SQL-based Mass Spectrometry Data Backend

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/MsqlBackend/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/MsqlBackend/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](http://codecov.io/github/RforMassSpectrometry/MsqlBackend/coverage.svg?branch=main)](http://codecov.io/github/RforMassSpectrometry/MsqlBackend?branch=main)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

This repository provides a *backend* for
[Spectra](https://github.com/RforMassSpectrometry/Spectra) objects that supports
storage of mass spectrometry (MS) data in an SQL database. The package provides
the functionality to create such databases from original (raw) MS data files (in
mzML, mzXML or netCDF format) and allows to extract the data in an efficient
way.

For more information see the package
[homepage](https://github.com/RforMassSpectrometry/MsqlBackend).

## Creating a database

By providing the connection to an SQL database, the `createMsqlBackendDatabase`
imports raw MS data from provided file names and stores it into the dedicated
database tables created during import. While `MsqlBackend` supports any type of
SQL database, it is currently optimized for MySQL/MariaDB databases.

## Using a *MsqlBackend* database

MS data in a *MsqlBackend* database can be accessed through the
[`Spectra`](https://github.com/RforMassSpectrometry/Spectra) package by using
the `MsqlBackend` MS backend. Assuming the variable `dbcon` represents a
(RDBI) database connection to a *MsqlBackend*, the data can be represented/used
with a `Spectra` object by:

```r
library(Spectra)
library(MsqlBackend)
sps <- Spectra(dbcon, source = MsqlBackend())
```
