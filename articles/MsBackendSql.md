# Storing Mass Spectrometry Data in SQL Databases

**Package**:
*[MsBackendSql](https://bioconductor.org/packages/3.23/MsBackendSql)*  
**Authors**: Johannes Rainer \[aut, cre\] (ORCID:
<https://orcid.org/0000-0002-6977-7147>), Chong Tang \[ctb\], Laurent
Gatto \[ctb\] (ORCID: <https://orcid.org/0000-0002-1520-2268>)  
**Compiled**: Thu Feb 5 10:12:25 2026

## Introduction

The *[Spectra](https://bioconductor.org/packages/3.23/Spectra)*
Bioconductor package provides a flexible and expandable infrastructure
for Mass Spectrometry (MS) data. The package supports interchangeable
use of different *backends* that provide additional file support or
different ways to store and represent MS data. The
*[MsBackendSql](https://bioconductor.org/packages/3.23/MsBackendSql)*
package provides backends to store data from whole MS experiments in SQL
databases. The data in such databases can be easily (and efficiently)
accessed using `Spectra` objects that use the `MsBackendSql` class as an
interface to the data in the database. Such `Spectra` objects have a
minimal memory footprint and hence allow analysis of very large data
sets even on computers with limited hardware capabilities. For certain
operations, the performance of this data representation is superior to
that of other low-memory (*on-disk*) data representations such as
`Spectra`’s `MsBackendMzR` backend. Finally, the `MsBackendSql` supports
also remote data access to e.g. a central database server hosting
several large MS data sets.

## Installation

The package can be installed with the `BiocManager` package. To install
`BiocManager` use `install.packages("BiocManager")` and, after that,
`BiocManager::install("MsBackendSql")` to install this package.

## Creating and using `MsBackendSql` SQL databases

`MsBackendSql` SQL databases can be created either by importing (raw) MS
data from MS data files using the
[`createMsBackendSqlDatabase()`](https://rformassspectrometry.github.io/MsBackendSql/reference/MsBackendSql.md)
or using the
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
function by providing in addition to the database connection also the
full MS data to import as a `DataFrame`. In the first example we use the
[`createMsBackendSqlDatabase()`](https://rformassspectrometry.github.io/MsBackendSql/reference/MsBackendSql.md)
function to import the full MS data from the provided MS data files into
an (empty) database. Below we first create an empty SQLite database (in
a temporary file) and use the
[`createMsBackendSqlDatabase()`](https://rformassspectrometry.github.io/MsBackendSql/reference/MsBackendSql.md)
function to create all necessary tables in that database and import the
MS data from two mzML files (provided throuth the
`r Biocpkg("MsDataHub")` package).

``` r

library(RSQLite)

dbfile <- tempfile()
con <- dbConnect(SQLite(), dbfile)

library(Spectra)
library(MsBackendSql)
fls <- c(MsDataHub::X20171016_POOL_POS_1_105.134.mzML(),
         MsDataHub::X20171016_POOL_POS_3_105.134.mzML())
createMsBackendSqlDatabase(con, fls)
dbDisconnect(con)
```

By default (with parameters `blob = TRUE` and
`peaksStorageMode = "blob2"`) the peaks data matrix of each spectrum is
stored as a *BLOB* data type into the database (one entry per spectrum).
This has advantages on the performance to extract the peaks data from
the database, but does not allow to filter individual peaks by their
*m/z* or intensity values directly in the database. As an alternative
(using `blob = FALSE`) it is also possible to store the individual *m/z*
and intensity values in separate columns of the database table. This
*long table format* results however in considerably larger databases
(with potentially poorer performance). Note also that the code and
backend is optimized for MySQL/MariaDB databases by taking advantage of
table partitioning and specialized table storage options. Any other SQL
database server is however also supported (also portable, self-contained
SQLite databases). In fact, performance for *MsBackendSql* databases
with peaks data stored as *BLOB* data type is similar for SQLite and
MySQL/MariaDB databases.

The *MsBackendSql* package provides two backends to interact with such
databases: the `MsBackendSql` class and the `MsBackendOfflineSql` class,
that inherits all properties and functions from the former, but does not
store the connection to the database within the object. The
`MsBackendOfflineSql` object thus supports parallel processing and
allows to save/load the object (e.g. using `save` and `saveRDS`). The
`MsBackendOfflineSql` might therefore be used as the preferred backend
to SQL databases for most applications.

To access the data in the database we create below a `Spectra` object
providing the database connection information in the constructor call
and specifying to use the `MsBackendOfflineSql` as *backend* (parameter
`source`). We stored the data to a SQLite database, thus we provide the
database name (SQLite database file name) and the SQLite DBI driver with
parameters `dbname` and `drv`. Which parameters are required to connect
to the database depends on the SQL database and the used driver. For a
MySQL/MariaDB database we would use the `MariaDB()` driver and would
have to provide the database name, user name, password as well as the
host name and port through which the database is accessible.

``` r

sps <- Spectra(dbname = dbfile, source = MsBackendOfflineSql(), drv = SQLite())
sps
```

    ## MSn data (Spectra) with 1862 spectra in a MsBackendOfflineSql backend:
    ##        msLevel precursorMz  polarity
    ##      <integer>   <numeric> <integer>
    ## 1            1          NA         1
    ## 2            1          NA         1
    ## 3            1          NA         1
    ## 4            1          NA         1
    ## 5            1          NA         1
    ## ...        ...         ...       ...
    ## 1858         1          NA         1
    ## 1859         1          NA         1
    ## 1860         1          NA         1
    ## 1861         1          NA         1
    ## 1862         1          NA         1
    ##  ... 35 more variables/columns.
    ##  Use  'spectraVariables' to list all of them.
    ## Database: /tmp/RtmpS7TBVj/file177c421a382

`Spectra` objects allow also to change the backend to any other backend
(extending `MsBackend`) using the
[`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
function. Below we use this function to first load all data into memory
by changing from the `MsBackendOfflineSql` to a `MsBackendMemory`.

``` r

sps_mem <- setBackend(sps, MsBackendMemory())
sps_mem
```

    ## MSn data (Spectra) with 1862 spectra in a MsBackendMemory backend:
    ##        msLevel     rtime scanIndex
    ##      <integer> <numeric> <integer>
    ## 1            1     0.280         1
    ## 2            1     0.559         2
    ## 3            1     0.838         3
    ## 4            1     1.117         4
    ## 5            1     1.396         5
    ## ...        ...       ...       ...
    ## 1858         1   258.636       927
    ## 1859         1   258.915       928
    ## 1860         1   259.194       929
    ## 1861         1   259.473       930
    ## 1862         1   259.752       931
    ##  ... 35 more variables/columns.
    ## Processing:
    ##  Switch backend from MsBackendOfflineSql to MsBackendMemory [Thu Feb  5 10:12:35 2026]

With this function it is also possible to change from any backend to a
`MsBackendOfflineSql` (or `MsBackendSql`) in which case a new database
is created and all data from the originating backend is stored in this
database. To change the backend to an `MsBackendOfflineSql` we need to
provide the connection information to the SQL database as additional
parameters. These parameters are the same that need to be passed to a
`dbConnect()` call to establish the connection to the database. These
parameters include the database driver (parameter `drv`), the database
name and eventually the user name, host etc (see `?dbConnect` for more
information). In the simple example below we store the data into a
SQLite database and thus only need to provide the database name, which
corresponds SQLite database file. In our example we store the data into
a temporary file. Optionally,
[`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
supports also the parameters `blob` and `peaksDataStorage` described
above for the
[`createMsBackendSqlDatabase()`](https://rformassspectrometry.github.io/MsBackendSql/reference/MsBackendSql.md)
function.

``` r

sps2 <- setBackend(sps_mem, MsBackendOfflineSql(), drv = SQLite(),
                   dbname = tempfile())
sps2
```

    ## MSn data (Spectra) with 1862 spectra in a MsBackendOfflineSql backend:
    ##        msLevel precursorMz  polarity
    ##      <integer>   <numeric> <integer>
    ## 1            1          NA         1
    ## 2            1          NA         1
    ## 3            1          NA         1
    ## 4            1          NA         1
    ## 5            1          NA         1
    ## ...        ...         ...       ...
    ## 1858         1          NA         1
    ## 1859         1          NA         1
    ## 1860         1          NA         1
    ## 1861         1          NA         1
    ## 1862         1          NA         1
    ##  ... 35 more variables/columns.
    ##  Use  'spectraVariables' to list all of them.
    ## Database: /tmp/RtmpS7TBVj/file177c3c78d238
    ## Processing:
    ##  Switch backend from MsBackendOfflineSql to MsBackendMemory [Thu Feb  5 10:12:35 2026]
    ##  Switch backend from MsBackendMemory to MsBackendOfflineSql [Thu Feb  5 10:12:35 2026]

Similar to any other `Spectra` object we can retrieve the available
*spectra variables* using the
[`spectraVariables()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
function.

``` r

spectraVariables(sps)
```

    ##  [1] "msLevel"                  "rtime"                   
    ##  [3] "acquisitionNum"           "scanIndex"               
    ##  [5] "dataStorage"              "dataOrigin"              
    ##  [7] "centroided"               "smoothed"                
    ##  [9] "polarity"                 "precScanNum"             
    ## [11] "precursorMz"              "precursorIntensity"      
    ## [13] "precursorCharge"          "collisionEnergy"         
    ## [15] "isolationWindowLowerMz"   "isolationWindowTargetMz" 
    ## [17] "isolationWindowUpperMz"   "peaksCount"              
    ## [19] "totIonCurrent"            "basePeakMZ"              
    ## [21] "basePeakIntensity"        "electronBeamEnergy"      
    ## [23] "ionisationEnergy"         "lowMZ"                   
    ## [25] "highMZ"                   "mergedScan"              
    ## [27] "mergedResultScanNum"      "mergedResultStartScanNum"
    ## [29] "mergedResultEndScanNum"   "injectionTime"           
    ## [31] "filterString"             "spectrumId"              
    ## [33] "ionMobilityDriftTime"     "scanWindowLowerLimit"    
    ## [35] "scanWindowUpperLimit"     "spectrum_id_"

The MS peak data can be accessed using either the
[`mz()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html),
[`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
or [`peaksData()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html)
functions. Below we extract the peaks matrix of the 5th spectrum and
display the first 6 rows.

``` r

peaksData(sps)[[5]] |>
head()
```

    ##            mz intensity
    ## [1,] 105.0347         0
    ## [2,] 105.0362       164
    ## [3,] 105.0376         0
    ## [4,] 105.0391         0
    ## [5,] 105.0405       328
    ## [6,] 105.0420         0

All data (peaks data or spectra variables) are **always** retrieved
on-the-fly from the database resulting thus in a minimal memory
footprint for the `Spectra` object.

``` r

print(object.size(sps), units = "KB")
```

    ## 116.3 Kb

The backend supports also adding additional spectra variables or
changing their values. Below we add 10 seconds to the retention time of
each spectrum.

``` r

sps$rtime <- sps$rtime + 10
```

Such operations do however **not** change the data in the database
(which is always considered read-only) but are cached locally within the
backend object (in memory). The size in memory of the object is thus
higher after changing that spectra variable.

``` r

print(object.size(sps), units = "KB")
```

    ## 131 Kb

Such `$<-` operations can also be used to *cache* spectra variables
(temporarily) in memory which can eventually improve performance. Below
we test the time it takes to extract the MS level from each spectrum
from the database, then cache the MS levels in memory using
`$msLevel <-` and test the timing to extract these cached variable.

``` r

system.time(msLevel(sps))
```

    ##    user  system elapsed 
    ##   0.010   0.000   0.009

``` r

sps$msLevel <- msLevel(sps)
system.time(msLevel(sps))
```

    ##    user  system elapsed 
    ##   0.005   0.000   0.004

We can also use the
[`reset()`](https://rdrr.io/pkg/Spectra/man/addProcessing.html) function
to *reset* the data to its original state (this will cause any local
spectra variables to be deleted and the backend to be initialized with
the original data in the database).

``` r

sps <- reset(sps)
```

## Performance considerations

### Database systems and data storage modes

The performance of storing and retrieving MS data from an `MsBackendSql`
respectively SQL database can also depend on the type of database used
as well as storage modes and the database layout used by `MsBackendSql`.

#### Database systems

Performance comparison have been made for small and large data sets
using different SQL database systems and *MsBackendSql* has been
optimized based on these results. For *MariaDB* database systems, for
example, the
[Aria](https://mariadb.com/docs/server/server-usage/storage-engines/aria)
storage engine is used by default as it has considerable advantages over
other MariaDB engines.

Performance of *MariaDB* and *SQLite* is comparable, even for very large
data sets/databases. See this [GitHub
issue](https://github.com/rformassspectrometry/MsBackendSql/issues/15)
for performance comparison between MariaDB and SQLite.

Performance evaluation of *SQLite* and *duckdb* are provided in this
[GitHub
issue](https://github.com/rformassspectrometry/MsBackendSql/issues/26).
*MsBackendSql* long format database layout (see next section for details
on available database layouts) with *duckdb* is clearly faster than with
*SQLite*. For the *blob2* database layout *SQLite* has advantages. Also,
extracting individual spectra variables or filtering by e.g. retention
time is slower for *duckdb*.

#### *MsBackendSql* database layouts/storage modes

*MsBackendSql* defines different database table layouts and hence ways
to store the MS data. The most intuitive way to store MS data would be
the *long* format (`peaksStorageMode = "long"`) which saves the *m/z*
and intensity values of each mass peak as a single row. While this would
allow to filter e.g. the peaks data by *m/z* and/or intensity values
already on the SQL level, it significantly increases the size of the
database. This is in particular true for *SQLite*-based databases. The
default storage mode (`peaksStorageMode = "blob2"`) stores the complete
peaks matrix (i.e. the two-column numerical matrix of *m/z* and
intensity values) of spectrum as one entity to the database. This entry
is stored as a binary data type (BLOB) in the database table (one row
per spectrum). This has a positive impact on the performance of the
database to extract peak data (which is much faster than from databases
with the *long* peaks storage mode). In addition, also the size (disk
space) of such databases are smaller. On the downside, these databases
will only be readable and usable with *MsBackendSql* or R-based tools.

For *MsBackendSql* in the *long* peaks storage mode it is suggested to
use *duckdb* as database backend.

### Performance comparison with other backends

The need to retrieve any spectra data on-the-fly from the database has
an impact on the performance of data access functions of `Spectra`
objects using `MsBackendSql`/`MsBackendOfflineSql` backends. To evaluate
this we compare below the performance of the `MsBackendSql` to other
`Spectra` backends, specifically, the `MsBackendMzR` which is the
default backend to read and represent raw MS data, and the
`MsBackendMemory` backend that keeps all MS data in memory (and is thus
not suggested for larger MS experiments). Similar to the `MsBackendMzR`,
also the `MsBackendSql` keeps only a limited amount of data in memory.
These *on-disk* backends need thus to retrieve spectra and MS peaks data
on-the-fly from either the original raw data files (in the case of the
`MsBackendMzR`) or from the SQL database (in the case of the
`MsBackendSql`). The in-memory backend `MsBackendMemory` is supposed to
provide the fastest data access since all data is kept in memory.

Below we thus create `Spectra` objects from the same data but using the
different backends.

``` r

con <- dbConnect(SQLite(), dbfile)
sps <- Spectra(con, source = MsBackendSql())
sps_mzr <- Spectra(fls, source = MsBackendMzR())
sps_im <- setBackend(sps_mzr, backend = MsBackendMemory())
```

At first we compare the memory footprint of the 3 backends.

``` r

print(object.size(sps), units = "KB")
```

    ## 114.7 Kb

``` r

print(object.size(sps_mzr), units = "KB")
```

    ## 401.1 Kb

``` r

print(object.size(sps_im), units = "KB")
```

    ## 54509 Kb

The `MsBackendSql` has the lowest memory footprint of all 3 backends
because it does not keep any data in memory. The `MsBackendMzR` keeps
all spectra variables, except the MS peaks data, in memory and has thus
a larger size. The `MsBackendMemory` keeps all data (including the MS
peaks data) in memory and has thus the largest size in memory.

Next we compare the performance to extract the MS level for each
spectrum from the 4 different `Spectra` objects.

``` r

library(microbenchmark)
microbenchmark(msLevel(sps),
               msLevel(sps_mzr),
               msLevel(sps_im))
```

    ## Unit: microseconds
    ##              expr      min       lq       mean    median       uq       max
    ##      msLevel(sps) 5092.230 5221.641 5527.46766 5344.7260 5706.677 10490.121
    ##  msLevel(sps_mzr)  479.305  508.479  539.91938  524.6735  575.213   764.076
    ##   msLevel(sps_im)   10.830   13.525   19.99937   20.4380   22.727    68.408
    ##  neval
    ##    100
    ##    100
    ##    100

Extracting MS levels is thus slowest for the `MsBackendSql`, which is
not surprising because both other backends keep this data in memory
while the `MsBackendSql` needs to retrieve it from the database.

We next compare the performance to access the full peaks data from each
`Spectra` object.

``` r

microbenchmark(peaksData(sps, BPPARAM = SerialParam()),
               peaksData(sps_mzr, BPPARAM = SerialParam()),
               peaksData(sps_im, BPPARAM = SerialParam()),
               times = 10)
```

    ## Unit: microseconds
    ##                                         expr        min         lq       mean
    ##      peaksData(sps, BPPARAM = SerialParam())  35421.545  40598.964 198973.941
    ##  peaksData(sps_mzr, BPPARAM = SerialParam()) 486077.553 487447.599 524937.162
    ##   peaksData(sps_im, BPPARAM = SerialParam())    418.812    532.904   1849.495
    ##       median         uq       max neval
    ##  200844.5155 355441.749 363044.14    10
    ##  490985.4995 500155.360 801863.45    10
    ##     666.9395    810.702  12765.87    10

As expected, the `MsBackendMemory` has the fasted access to the full
peaks data. The `MsBackendSql` outperforms however the `MsBackendMzR`
providing faster access to the m/z and intensity values.

Performance can be improved for the `MsBackendMzR` using parallel
processing. Note that the `MsBackendSql` does **not support** parallel
processing and thus parallel processing is (silently) disabled in
functions such as
[`peaksData()`](https://rdrr.io/pkg/ProtGenerics/man/peaksData.html).

``` r

m2 <- MulticoreParam(2)
microbenchmark(peaksData(sps, BPPARAM = m2),
               peaksData(sps_mzr, BPPARAM = m2),
               peaksData(sps_im, BPPARAM = m2),
               times = 10)
```

    ## Unit: microseconds
    ##                              expr        min        lq        mean     median
    ##      peaksData(sps, BPPARAM = m2)  35999.834  48735.94 100612.6862  60387.819
    ##  peaksData(sps_mzr, BPPARAM = m2) 420067.129 433398.50 689278.4337 476967.549
    ##   peaksData(sps_im, BPPARAM = m2)    552.701    809.21    887.4869    837.262
    ##          uq         max neval
    ##   68584.496  486145.978    10
    ##  979909.664 1336642.221    10
    ##     954.171    1241.797    10

We next compare the performance of subsetting operations.

``` r

microbenchmark(filterRt(sps, rt = c(50, 100)),
               filterRt(sps_mzr, rt = c(50, 100)),
               filterRt(sps_im, rt = c(50, 100)))
```

    ## Unit: microseconds
    ##                                expr      min        lq      mean   median
    ##      filterRt(sps, rt = c(50, 100)) 1831.316 1850.3175 2037.5777 1875.705
    ##  filterRt(sps_mzr, rt = c(50, 100)) 1348.105 1392.2725 1514.9231 1421.597
    ##   filterRt(sps_im, rt = c(50, 100))  399.806  426.9065  451.3687  445.361
    ##         uq       max neval
    ##  1903.2965 14769.796   100
    ##  1448.4125 10039.922   100
    ##   458.3155  1077.419   100

The two *on-disk* backends `MsBackendSql` and `MsBackendMzR` show a
comparable performance for this operation. This filtering does involves
access to a spectra variables (the retention time in this case) which,
for the `MsBackendSql` needs first to be retrieved from the backend. The
`MsBackendSql` backend allows however also to *cache* spectra variables
(i.e. they are stored within the `MsBackendSql` object). Any access to
such cached spectra variables can eventually be faster because no
dedicated SQL query is needed.

To evaluate the performance of a *pure* subsetting operation we first
define the indices of 10 random spectra and subset the `Spectra` objects
to these.

``` r

idx <- sample(seq_along(sps), 10)
microbenchmark(sps[idx],
               sps_mzr[idx],
               sps_im[idx])
```

    ## Unit: microseconds
    ##          expr     min       lq     mean  median       uq      max neval
    ##      sps[idx] 136.054 142.6910 153.5289 153.281 160.2840  205.163   100
    ##  sps_mzr[idx] 664.550 683.9360 709.1889 692.597 704.7050 2075.402   100
    ##   sps_im[idx] 241.060 250.2865 260.6238 259.429 266.2515  348.620   100

Here the `MsBackendSql` outperforms the other backends because it does
not keep any data in memory and hence does not need to subset these. The
two other backends need to subset the data they keep in memory which is
in both cases a data frame with either a reduced set of spectra
variables or the full MS data.

At last we compare also the extraction of the peaks data from the such
subset `Spectra` objects.

``` r

sps_10 <- sps[idx]
sps_mzr_10 <- sps_mzr[idx]
sps_im_10 <- sps_im[idx]

microbenchmark(peaksData(sps_10),
               peaksData(sps_mzr_10),
               peaksData(sps_im_10),
               times = 10)
```

    ## Unit: microseconds
    ##                   expr       min        lq      mean    median        uq
    ##      peaksData(sps_10)  1982.218  2129.522  2674.657  2381.002  3302.902
    ##  peaksData(sps_mzr_10) 57582.971 60943.230 62377.825 61947.523 63673.113
    ##   peaksData(sps_im_10)   405.737   459.908   584.789   631.989   675.229
    ##        max neval
    ##   3497.324    10
    ##  68307.109    10
    ##    752.464    10

The `MsBackendSql` outperforms the `MsBackendMzR` while, not
unexpectedly, the `MsBackendMemory` provides fasted access.

### Considerations for database systems/servers

The backends from the *MsBackendSql* package use standard SQL calls to
retrieve MS data from the database and hence any SQL database system
(for which an R package is available) is supported. SQLite-based
databases would represent the easiest and most user friendly solution
since no database server administration and user management is required.
Indeed, performance of SQLite is very high, even for very large data
sets. Server-based databases on the other hand have the advantage to
enable a centralized storage and control of MS data (inclusive user
management etc). Also, such server systems would also allow data set or
server-specific configurations to improve performance.

A comparison between a SQLite-based with a MariaDB-based *MsBackendSql*
database for a large data set comprising over 8,000 samples and over
15,000,000 spectra is available
[here](https://github.com/rformassspectrometry/MsBackendSql/issues/15).
In brief, performance to extract data was comparable and for individual
spectra variables even faster for the SQLite database. Only when more
complex SQL queries were involved (combining several primary keys or
data fields) the more advanced MariaDB database outperformed SQLite.

## Other properties of the `MsBackendSql`

The `MsBackendSql` backend does not support parallel processing since
the database connection can not be shared across the different
(parallel) processes. Thus, all methods on `Spectra` objects that use a
`MsBackendSql` will automatically (and silently) disable parallel
processing even if a dedicated parallel processing setup was passed
along with the `BPPARAM` method.

Some functions on `Spectra` objects require to load the MS peak data
(i.e., m/z and intensity values) into memory. For very large data sets
(or computers with limited hardware resources) such function calls can
cause out-of-memory errors. One example is the
[`lengths()`](https://rdrr.io/r/base/lengths.html) function that
determines the number of peaks per spectrum by loading the peak matrix
first into memory. Such functions should ideally be called using the
`peaksapply()` function with parameter `chunkSize` (e.g.,
`peaksapply(sps, lengths, chunkSize = 5000L)`). Instead of processing
the full data set, the data will be first split into chunks of size
`chunkSize` that are stepwise processed. Hence, only data from
`chunkSize` spectra is loaded into memory in one iteration.

## Summary

The `MsBackendSql` provides an MS data representations and storage mode
with a minimal memory footprint (in R) that is still comparably
efficient for standard processing and subsetting operations. This
backend is specifically useful for very large MS data sets, that could
even be hosted on remote (MySQL/MariaDB) servers. A potential use case
for this backend could thus be to set up a central storage place for MS
experiments with data analysts connecting remotely to this server to
perform initial data exploration and filtering. After subsetting to a
smaller data set of interest, users could then retrieve/download this
data by changing the backend to e.g. a `MsBackendMemory`, which would
result in a *download* of the full data to the user computer’s memory.

## Session information

``` r

sessionInfo()
```

    ## R Under development (unstable) (2026-02-01 r89366)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] microbenchmark_1.5.0 MsDataHub_1.11.0     RSQLite_2.4.5       
    ##  [4] MsBackendSql_1.11.3  Spectra_1.21.1       BiocParallel_1.45.0 
    ##  [7] S4Vectors_0.49.0     BiocGenerics_0.57.0  generics_0.1.4      
    ## [10] BiocStyle_2.39.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1       dplyr_1.2.0            blob_1.3.0            
    ##  [4] filelock_1.0.3         Biostrings_2.79.4      fastmap_1.2.0         
    ##  [7] BiocFileCache_3.1.0    digest_0.6.39          lifecycle_1.0.5       
    ## [10] cluster_2.1.8.2        ProtGenerics_1.43.0    KEGGREST_1.51.1       
    ## [13] magrittr_2.0.4         compiler_4.6.0         rlang_1.1.7           
    ## [16] sass_0.4.10            progress_1.2.3         tools_4.6.0           
    ## [19] yaml_2.3.12            data.table_1.18.2.1    knitr_1.51            
    ## [22] prettyunits_1.2.0      htmlwidgets_1.6.4      bit_4.6.0             
    ## [25] curl_7.0.0             withr_3.0.2            purrr_1.2.1           
    ## [28] desc_1.4.3             ExperimentHub_3.1.0    MASS_7.3-65           
    ## [31] cli_3.6.5              mzR_2.45.0             rmarkdown_2.30        
    ## [34] crayon_1.5.3           ragg_1.5.0             otel_0.2.0            
    ## [37] httr_1.4.7             ncdf4_1.24             DBI_1.2.3             
    ## [40] cachem_1.1.0           parallel_4.6.0         AnnotationDbi_1.73.0  
    ## [43] BiocManager_1.30.27    XVector_0.51.0         vctrs_0.7.1           
    ## [46] jsonlite_2.0.0         bookdown_0.46          IRanges_2.45.0        
    ## [49] hms_1.1.4              bit64_4.6.0-1          clue_0.3-66           
    ## [52] systemfonts_1.3.1      jquerylib_0.1.4        glue_1.8.0            
    ## [55] pkgdown_2.2.0.9000     codetools_0.2-20       stringi_1.8.7         
    ## [58] BiocVersion_3.23.1     tibble_3.3.1           pillar_1.11.1         
    ## [61] rappdirs_0.3.4         htmltools_0.5.9        Seqinfo_1.1.0         
    ## [64] R6_2.6.1               dbplyr_2.5.1           httr2_1.2.2           
    ## [67] textshaping_1.0.4      evaluate_1.0.5         Biobase_2.71.0        
    ## [70] AnnotationHub_4.1.0    png_0.1-8              memoise_2.0.1         
    ## [73] bslib_0.10.0           MetaboCoreUtils_1.19.1 Rcpp_1.1.1            
    ## [76] fastmatch_1.1-8        xfun_0.56              MsCoreUtils_1.23.2    
    ## [79] fs_1.6.6               pkgconfig_2.0.3
