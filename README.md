# SQL-based Mass Spectrometry Data Backend

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/MsBackendSql/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/MsBackendSql/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov](https://codecov.io/gh/rformassspectrometry/MsBackendSql/branch/devel/graph/badge.svg?token=9qCNOICdYv)](https://codecov.io/gh/rformassspectrometry/MsBackendSql)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)
[![years in bioc](http://bioconductor.org/shields/years-in-bioc/MsBackendSql.svg)](https://bioconductor.org/packages/release/bioc/html/MsBackendSql.html)
[![Ranking by downloads](http://bioconductor.org/shields/downloads/release/MsBackendSql.svg)](https://bioconductor.org/packages/stats/bioc/MsBackendSql/)
[![build release](http://bioconductor.org/shields/build/release/bioc/MsBackendSql.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MsBackendSql/)
[![build devel](http://bioconductor.org/shields/build/devel/bioc/MsBackendSql.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/MsBackendSql/)

This repository provides a *backend* for
[Spectra](https://github.com/RforMassSpectrometry/Spectra) objects that supports
storage of mass spectrometry (MS) data in an SQL database. The package provides
the functionality to create such databases from original (raw) MS data files (in
mzML, mzXML or netCDF format) and allows to extract the data in an efficient
way.

For more information see the package
[homepage](https://github.com/RforMassSpectrometry/MsBackendSql).

## Creating a database

By providing the connection to an SQL database, the `createMsBackendSqlDatabase`
imports raw MS data from provided file names and stores it into the dedicated
database tables created during import. While `MsBackendSql` supports any type of
SQL database, it is currently optimized for MySQL/MariaDB databases.

## Using a *MsBackendSql* database

MS data in a *MsBackendSql* database can be accessed through the
[`Spectra`](https://github.com/RforMassSpectrometry/Spectra) package by using
the `MsBackendSql` MS backend. Assuming the variable `dbcon` represents a
(RDBI) database connection to a *MsBackendSql*, the data can be represented/used
with a `Spectra` object by:

```r
library(Spectra)
library(MsBackendSql)
sps <- Spectra(dbcon, source = MsBackendSql())
```
