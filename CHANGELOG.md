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
