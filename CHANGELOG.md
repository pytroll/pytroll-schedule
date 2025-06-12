## Version 0.7.3 (2025/04/30)

### Issues Closed

* [Issue 99](https://github.com/pytroll/pytroll-schedule/issues/99) - Build sometimes allocates 24 GB of RAM ([PR 100](https://github.com/pytroll/pytroll-schedule/pull/100) by [@mraspaud](https://github.com/mraspaud))
* [Issue 96](https://github.com/pytroll/pytroll-schedule/issues/96) - Add fengyun3f mersi-3
* [Issue 65](https://github.com/pytroll/pytroll-schedule/issues/65) - Mapping in collections is depricated since 3.3 and removed in 3.10

In this release 3 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 100](https://github.com/pytroll/pytroll-schedule/pull/100) - Fix test running long ([99](https://github.com/pytroll/pytroll-schedule/issues/99))

In this release 1 pull request was closed.


###############################################################################
## Version 0.7.2 (2025/01/08)

### Issues Closed

* [Issue 90](https://github.com/pytroll/pytroll-schedule/issues/90) - Should be possible to store swath outlines as shapefiles ([PR 93](https://github.com/pytroll/pytroll-schedule/pull/93) by [@adybbroe](https://github.com/adybbroe))
* [Issue 88](https://github.com/pytroll/pytroll-schedule/issues/88) - Change the avoid command line options to take a list ([PR 89](https://github.com/pytroll/pytroll-schedule/pull/89) by [@TAlonglong](https://github.com/TAlonglong))
* [Issue 69](https://github.com/pytroll/pytroll-schedule/issues/69) - Add slstr as instrument ([PR 70](https://github.com/pytroll/pytroll-schedule/pull/70) by [@TAlonglong](https://github.com/TAlonglong))

In this release 3 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 92](https://github.com/pytroll/pytroll-schedule/pull/92) - Bugfix - slstr not yet properly supported
* [PR 86](https://github.com/pytroll/pytroll-schedule/pull/86) - Add a test for runing from start to finish
* [PR 81](https://github.com/pytroll/pytroll-schedule/pull/81) - Fix ftp downloading for older python versions

#### Features added

* [PR 95](https://github.com/pytroll/pytroll-schedule/pull/95) - Fix a few minor issues raised in PR93 which were accidentally left out
* [PR 93](https://github.com/pytroll/pytroll-schedule/pull/93) - Create instrument swath outlines as shapefiles ([90](https://github.com/pytroll/pytroll-schedule/issues/90))
* [PR 89](https://github.com/pytroll/pytroll-schedule/pull/89) - Issue 88 command line option avoid as list ([88](https://github.com/pytroll/pytroll-schedule/issues/88))
* [PR 87](https://github.com/pytroll/pytroll-schedule/pull/87) - Add noaa21
* [PR 86](https://github.com/pytroll/pytroll-schedule/pull/86) - Add a test for runing from start to finish
* [PR 70](https://github.com/pytroll/pytroll-schedule/pull/70) - Issue 69 add slstr ([69](https://github.com/pytroll/pytroll-schedule/issues/69))

In this release 9 pull requests were closed.


## Version 0.7.1 (2024/02/16)


### Pull Requests Merged

#### Features added

* [PR 85](https://github.com/pytroll/pytroll-schedule/pull/85) - Update CI Python versions
* [PR 84](https://github.com/pytroll/pytroll-schedule/pull/84) - Update versioneer

In this release 2 pull requests were closed.

## Version 0.7.0 (2023/06/21)

### Issues Closed

* [Issue 74](https://github.com/pytroll/pytroll-schedule/issues/74) - numpy 1.24.0 np.bool removed needs cleanup ([PR 72](https://github.com/pytroll/pytroll-schedule/pull/72) by [@mraspaud](https://github.com/mraspaud))
* [Issue 68](https://github.com/pytroll/pytroll-schedule/issues/68) - Test failure with PyResample 1.23 ([PR 72](https://github.com/pytroll/pytroll-schedule/pull/72) by [@mraspaud](https://github.com/mraspaud))

In this release 2 issues were closed.

### Pull Requests Merged

#### Features added

* [PR 72](https://github.com/pytroll/pytroll-schedule/pull/72) - Start refactoring pytroll schedule ([74](https://github.com/pytroll/pytroll-schedule/issues/74), [68](https://github.com/pytroll/pytroll-schedule/issues/68))
* [PR 71](https://github.com/pytroll/pytroll-schedule/pull/71) - Create dependabot.yml

In this release 2 pull requests were closed.


## Version 0.6.0 (2021/12/09)

### Issues Closed

* [Issue 62](https://github.com/pytroll/pytroll-schedule/issues/62) - Remove remnants of Python 2 support ([PR 67](https://github.com/pytroll/pytroll-schedule/pull/67) by [@pnuu](https://github.com/pnuu))
* [Issue 60](https://github.com/pytroll/pytroll-schedule/issues/60) - Deprecated import of Mapping
* [Issue 59](https://github.com/pytroll/pytroll-schedule/issues/59) - Failures in Schedule tests ([PR 61](https://github.com/pytroll/pytroll-schedule/pull/61) by [@pnuu](https://github.com/pnuu))
* [Issue 54](https://github.com/pytroll/pytroll-schedule/issues/54) - Deprecated use of abstract base classes ([PR 57](https://github.com/pytroll/pytroll-schedule/pull/57) by [@pnuu](https://github.com/pnuu))
* [Issue 53](https://github.com/pytroll/pytroll-schedule/issues/53) - The unittests are not run automatically ([PR 55](https://github.com/pytroll/pytroll-schedule/pull/55) by [@pnuu](https://github.com/pnuu))
* [Issue 52](https://github.com/pytroll/pytroll-schedule/issues/52) - Boundary calculations are broken ([PR 56](https://github.com/pytroll/pytroll-schedule/pull/56) by [@pnuu](https://github.com/pnuu))
* [Issue 49](https://github.com/pytroll/pytroll-schedule/issues/49) - Three unit tests failed.

In this release 7 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 61](https://github.com/pytroll/pytroll-schedule/pull/61) - Allow `mersi-2` as instrument name ([59](https://github.com/pytroll/pytroll-schedule/issues/59))
* [PR 56](https://github.com/pytroll/pytroll-schedule/pull/56) - Remove a bug introduced in PR38 ([52](https://github.com/pytroll/pytroll-schedule/issues/52))
* [PR 51](https://github.com/pytroll/pytroll-schedule/pull/51) - Remove some redundant code and fix a failed unit test.
* [PR 45](https://github.com/pytroll/pytroll-schedule/pull/45) - Use recent ssl protocol for older python versions
* [PR 38](https://github.com/pytroll/pytroll-schedule/pull/38) - Fix S3 olci scan duration

#### Features added

* [PR 67](https://github.com/pytroll/pytroll-schedule/pull/67) - Refactor remove legacy code support ([62](https://github.com/pytroll/pytroll-schedule/issues/62))
* [PR 66](https://github.com/pytroll/pytroll-schedule/pull/66) - Change tested Python versions to 3.8, 3.9 and 3.10
* [PR 64](https://github.com/pytroll/pytroll-schedule/pull/64) - Use safe loading for YAML config file
* [PR 61](https://github.com/pytroll/pytroll-schedule/pull/61) - Allow `mersi-2` as instrument name ([59](https://github.com/pytroll/pytroll-schedule/issues/59))
* [PR 58](https://github.com/pytroll/pytroll-schedule/pull/58) - Fix a test failure on Python 3.7
* [PR 57](https://github.com/pytroll/pytroll-schedule/pull/57) - Fix an import raising deprecation warning ([54](https://github.com/pytroll/pytroll-schedule/issues/54))
* [PR 55](https://github.com/pytroll/pytroll-schedule/pull/55) - Add GitHub actions to run unittests ([53](https://github.com/pytroll/pytroll-schedule/issues/53))
* [PR 50](https://github.com/pytroll/pytroll-schedule/pull/50) - Add a southern hemisphere pass test.
* [PR 46](https://github.com/pytroll/pytroll-schedule/pull/46) - Give the option to plot multiple polygons
* [PR 45](https://github.com/pytroll/pytroll-schedule/pull/45) - Use recent ssl protocol for older python versions
* [PR 44](https://github.com/pytroll/pytroll-schedule/pull/44) - Make plot filename more complete, including the instrument name
* [PR 42](https://github.com/pytroll/pytroll-schedule/pull/42) - Make it possible to tell cartopy to use offline shapefiles
* [PR 41](https://github.com/pytroll/pytroll-schedule/pull/41) - Fix nasa ftp retrieval
* [PR 38](https://github.com/pytroll/pytroll-schedule/pull/38) - Fix S3 olci scan duration

In this release 19 pull requests were closed.


## Version 0.5.2 (2019/03/19)


### Pull Requests Merged

#### Bugs fixed

* [PR 36](https://github.com/pytroll/pytroll-schedule/pull/36) - Add xarray to conda dependencies
* [PR 35](https://github.com/pytroll/pytroll-schedule/pull/35) - Bugfix - when a set of sensors are provided choose avhrr if it is one of them
* [PR 33](https://github.com/pytroll/pytroll-schedule/pull/33) - Bugfix avhrr naming

In this release 3 pull requests were closed.

## Version 0.5.1 (2019/01/08)

### Issues Closed

* [Issue 27](https://github.com/pytroll/pytroll-schedule/issues/27) - Drawing the ascat outline
* [Issue 25](https://github.com/pytroll/pytroll-schedule/issues/25) - New version slower in generating the schedule ([PR 26](https://github.com/pytroll/pytroll-schedule/pull/26))

In this release 2 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 26](https://github.com/pytroll/pytroll-schedule/pull/26) - Speed up schedule generation ([25](https://github.com/pytroll/pytroll-schedule/issues/25))
* [PR 24](https://github.com/pytroll/pytroll-schedule/pull/24) - Bugfix schedule generation

#### Features added

* [PR 29](https://github.com/pytroll/pytroll-schedule/pull/29) - Move save_fig import in under the save function
* [PR 28](https://github.com/pytroll/pytroll-schedule/pull/28) - Restructure the get_next_passes putting some of the highly nested cod…
* [PR 26](https://github.com/pytroll/pytroll-schedule/pull/26) - Speed up schedule generation ([25](https://github.com/pytroll/pytroll-schedule/issues/25))
* [PR 23](https://github.com/pytroll/pytroll-schedule/pull/23) - Versioneer

In this release 6 pull requests were closed.

## Version 0.5.0 (2018/11/25)

### Issues Closed

* [Issue 17](https://github.com/pytroll/pytroll-schedule/issues/17) - Plotting facility uses basemap ([PR 18](https://github.com/pytroll/pytroll-schedule/pull/18))

In this release 1 issue was closed.
