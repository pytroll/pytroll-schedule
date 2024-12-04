Changelog
=========


v0.4.0 (2018-11-08)
-------------------

Fix
~~~
- Bugfix: Resolve import errors - a few things have been moved to
  pyresample. [Adam.Dybbroe]
- Bugfix: setting the instrument dependent scan duration. [Adam.Dybbroe]
- Bugfix: default instrument geometry function. [Adam.Dybbroe]
- Bugfix: Instrument for Suomi NPP was not known. [Adam.Dybbroe]
- Multi-proc only if stations>1 ; eval Aqua-dump with default URL.
  [Alexander Maul]

Other
~~~~~
- Update changelog. [Martin Raspaud]
- Bump version: 0.3.3 → 0.4.0. [Martin Raspaud]
- Merge pull request #16 from pytroll/debug-satpass. [Martin Raspaud]

  Debug satpass
- Make sure indecies are integers. [Adam.Dybbroe]
- Fix boundary import from pyresample. [Martin Raspaud]
- Merge branch 'debug-satpass' of github.com:pytroll/pytroll-schedule
  into debug-satpass. [Adam.Dybbroe]
- Use pyresample's boundary module. [Martin Raspaud]
- Add unittest module for satpass testing. [Adam.Dybbroe]
- Add unittests for swath boundary and swath coverage. [Adam.Dybbroe]
- Remove "v" from versioning. [Adam.Dybbroe]
- Bugfix bottom and top lons,lats affecting at least VIIRS.
  [Adam.Dybbroe]
- Remove log info instrument Signed-off-by: Adam.Dybbroe
  <adam.dybbroe@smhi.se> [Adam.Dybbroe]
- Merge branch 'master' into debug-satpass. [Adam.Dybbroe]
- Merge pull request #11 from pytroll/feature-python3-support. [Panu
  Lahtinen]

  Support Python 3
- Encapsulate direct call of the files inside main() function. [Panu
  Lahtinen]
- Remove leftover rebase marker. [Panu Lahtinen]
- Merge branch 'feature-python3-support' of https://github.com/mraspaud
  /pytroll-schedule into feature-python3-support. [Panu Lahtinen]
- Merge branch 'feature-python3-support' of https://github.com/mraspaud
  /pytroll-schedule into feature-python3-support. [Panu Lahtinen]
- Fix urlparse imports. [Panu Lahtinen]
- Remove unused import. [Panu Lahtinen]
- Support Python 3. [Panu Lahtinen]
- Fix urlparse imports. [Panu Lahtinen]
- Remove unused import. [Panu Lahtinen]
- Support Python 3. [Panu Lahtinen]
- Support Python 3. [Panu Lahtinen]
- Fix urlparse imports. [Panu Lahtinen]
- Remove unused import. [Panu Lahtinen]
- Support Python 3. [Panu Lahtinen]
- Bugfix checking instrument name. [Adam.Dybbroe]
- Pep8. [Adam.Dybbroe]
- Fix sec_scan_duration for avhrr and ascat. [Adam.Dybbroe]
- Use variable scan_step instead of frequency. [Adam.Dybbroe]
- Check if instrument is a list. [Adam.Dybbroe]
- A list of instruments is not allowed. Set to avhrr for the time being.
  [Adam.Dybbroe]
- Bugfix. [Adam.Dybbroe]
- Add more debug info... [Adam.Dybbroe]
- Add debug info. [Adam.Dybbroe]
- Merge pull request #15 from pytroll/debug-schedule-page-generation.
  [Adam Dybbroe]

  Debug schedule page generation
- Bugfix viirs and getting the sides right. [Adam.Dybbroe]
- Add more debug info. [Adam.Dybbroe]
- Fix exception message and add debug printouts. [Adam.Dybbroe]
- Merge pull request #14 from pytroll/bugfix-snpp. [Adam Dybbroe]

  Bugfix: Instrument for Suomi NPP was not known.
- Merge pull request #12 from pytroll/feature-ascat. [Adam Dybbroe]

  Add support for ASCAT scan geometry
- Bugfix, set the correct instrument when creating the Pass object.
  [Adam.Dybbroe]
- Fix for Python3 and add number of FOVs for mhs and amsua.
  [Adam.Dybbroe]
- Generalize to other instruments different from AVHRR. [Adam.Dybbroe]
- Add support for ASCAT scan geometry. [Adam.Dybbroe]
- Merge pull request #10 from pytroll/develop. [Martin Raspaud]

  Get rid of develop
- Merge pull request #9 from pytroll/feature-oo. [Martin Raspaud]

  [WIP] Factorize code using OOP
- Adapt legacy cfg reader to use the new classes. [Alexander Maul]
- Make more code working in Py2.7 & Py3.4+. [Alexander Maul]
- Update documentation for config file (plain->yaml). [Alexander Maul]
- Remove surplus module function and parameter, both now in classes.
  [Alexander Maul]
- Use AQUA/TERRA dump URL from yaml-cfg, add use_2to3 to setup.
  [Alexander Maul]
- Merge branch 'feature-oo' of https://github.com/pytroll/pytroll-
  schedule.git into feature-oo. [Alexander Maul]
- Ammend yaml-config/reader to set sat-scores per station (optional).
  [Alexander Maul]
- Ammend yaml-config/reader to set sat-scores per station (optional).
  [Alexander Maul]
- Change print() statements to proper logging in exception handling.
  [Alexander Maul]
- Change combined_stations() for oo-rewrite, correct some typos for py3.
  [Alexander Maul]
- Fix circular import. [Martin Raspaud]
- Add support for olci instrument. [Martin Raspaud]
- Add python 3 support. [Martin Raspaud]
- Restore string as satellite for pass instatiation. [Martin Raspaud]
- Change default mail sender. [Martin Raspaud]
- Create some classes and use them. [Martin Raspaud]
- Let open generate an error if the config file is missing. [Martin
  Raspaud]
- Add more option to the drawing feature. [Martin Raspaud]
- Bugfix for metop-a/b and noaa-20. [Adam.Dybbroe]
- Correct NOAA-20 naming for TLEs, in case JPSS-1 is used instead of
  "NOAA 20" [Adam.Dybbroe]
- Bugfix - exclude satellites. [Adam.Dybbroe]
- Allow to exclude satellite platforms via the command line.
  [Adam.Dybbroe]
- Fix for NOAA 20 / JPSS-1 naming. [Adam.Dybbroe]
- Handle situations more gracefully with satellites for which we have no
  TLEs. And prepare for NOAA-20. [Adam.Dybbroe]
- Fix center_id for combined schedules. [Alexander Maul]
- Merge pull request #8 from pytroll/yaml-config. [Alex Maul]

  Yaml config with requested changes.
- Merge branch 'develop' into yaml-config. [Alex Maul]
- Make center_id a configuration item. [Martin Raspaud]
- Update usage in docs. [Martin Raspaud]
- Merge branch 'develop' of github.com:mraspaud/pytroll-schedule into
  develop. [Martin Raspaud]
- Allow passing proj info to the mapper. [Martin Raspaud]
- Update with changes according to PR reviewer's comments. [Alexander
  Maul]

  Change key "center" to "center_id" yaml-config,
  move lines evaluating new keys up in trollsched.utils.
- Move center-ID to yaml-config. [Alexander Maul]
- Add "minimum pass duration" and "dump host URL" to yaml config file.
  [Alexander Maul]
- Merge pull request #7 from pytroll/develop. [Alex Maul]

  Include last commits into branch yaml-config.
- Fix filename filter for Aqua/Terra transponder-off info files.
  [Alexander Maul]
- Merge branch 'develop' into develop. [Alex Maul]
- Merge pull request #4 from pytroll/develop. [Alex Maul]

  Update fork
- Increase estimate for the size of the combined passes' graph.
  [Alexander Maul]
- Move code collecting labels per tree-node into local function.
  [Alexander Maul]
- Remove obsolete "to-do" comment. [Alexander Maul]
- Improve __eq__ in SimplePass. [Alexander Maul]

  Now, passes are compared by name and orbit number if both are instances
  of Pass.
  Otherwise the old comparision is used.

- Add missing import. [Alexander Maul]
- Merge pull request #2 from pytroll/develop. [Alex Maul]

  update develop-branch at alexmaul from pytroll
- Update docs/gitignore. [Alexander Maul]
- Add information on Direct Broadcast Downlink Scheduling restrictions
  for EOS and Metop. [Adam.Dybbroe]
- Add YAML reader for configuration file. [Alexander Maul]

  Both the old config-format and yaml should be supported.
  Also add a template for a yaml configuration file.
- Yaml-config reader: start work. [alexmaul]


v0.3.3 (2017-09-20)
-------------------
- Update changelog. [Martin Raspaud]
- Bump version: 0.3.2 → 0.3.3. [Martin Raspaud]
- Remove support for 2.6. [Martin Raspaud]
- Add/change function to load/include info about TERRA dumps, in
  addition to AQUA dumps. [alexmaul]


v0.3.2 (2017-08-18)
-------------------

Fix
~~~
- Bugfix: Number of scans should be an integer. [Adam.Dybbroe]
- Bugfix: Correct path to output png plots. [Adam.Dybbroe]

Other
~~~~~
- Update changelog. [Martin Raspaud]
- Bump version: 0.3.1 → 0.3.2. [Martin Raspaud]
- Merge remote-tracking branch 'origin/generate-schedule-pages' into
  develop. [Martin Raspaud]
- Use TLES environment varaible. [Adam.Dybbroe]
- Clean up and make sat names look nice on plots. [Adam.Dybbroe]
- Add a schedule generator runner based on posttroll messaging.
  [Adam.Dybbroe]
- Fix int division for py3. [Martin Raspaud]
- Fix unit test test_normalize() [Panu Lahtinen]
- Add IASI scan angle to instrument swath boundary calculations. [Panu
  Lahtinen]
- Merge pull request #4 from alexmaul/develop. [Martin Raspaud]

  Develop
- First draft on scheduler documentation. [Alexander Maul]
- Fix typo in send_file() [Alexander Maul]
- Collected all styles in head/style and fixed font-sizes. [Alexander
  Maul]
- Always save xml-file in request-mode. [Alexander Maul]

  Even if report-mode is set in command-line arguments an xml-file in
  request-mode is created.
  Also the combined-request files are transfered with FTP, which is
  encapsuled in a new function.

- Use div with css-positioning instead of sturdy tables. [Alexander
  Maul]
- Create XSL for display of the aquisition-report in a browser.
  [Alexander Maul]
- Version-compatible dictionary building list-comprehension. [Alexander
  Maul]
- Fix missing setter. [Alexander Maul]
- Merge pull request #1 from pytroll/develop. [Alex Maul]

  Sync Develop to fork
- Corrected last typos. [Alexander Maul]
- Merge branch 'multiple_stations' into develop. [Alexander Maul]

  # Conflicts:
  #	trollsched/schedule.py

- All test work fine, prepare for PR. [Alexander Maul]
- Passes handling and other re-works finished. [Alexander Maul]
- Intermediate commit. [Alexander Maul]
- Intermediate commit. combining three stations works. [Alexander Maul]
- Last quirks removed ... intense testing follows. [Alexander Maul]
- Intermediate commit to backup work. [Alexander Maul]
- Intermediate commit. [Alexander Maul]

  it's working. generated path is indeed more optimal.
- Combining weighted trees works so far. not ready for pull-req.
  [Alexander Maul]
- Schedule.py : remove unneccessary writing to temp-file. [Alexander
  Maul]

  combine.py : just saving some development, nothing ready yet.
- Most of the changes from workshop. Starts branch in own fork.
  [Alexander Maul]
- Intergrate changes from WS. [Alexander Maul]
- Index on develop: 40a9016 Add avoid_list feature. [Alexander Maul]


v0.3.1 (2016-12-08)
-------------------
- Update changelog. [Martin Raspaud]
- Bump version: 0.3.0 → 0.3.1. [Martin Raspaud]
- Bugfix in spherical. [Martin Raspaud]


v0.3.0 (2016-10-27)
-------------------

Changes
~~~~~~~
- Allow to pass a parameter to modpi. [Martin Raspaud]

Fix
~~~
- Bugfix: Don't duplicate points in decimating boundaries. [Martin
  Raspaud]
- Bugfix: rec flag cannot be bool for xml. [Martin Raspaud]
- Bugfix: save rec status, not image link... [Martin Raspaud]

Other
~~~~~
- Update changelog. [Martin Raspaud]
- Bump version: 0.2.2 → 0.3.0. [Martin Raspaud]
- Add bump and changelog configs. [Martin Raspaud]
- Simplify version management. [Martin Raspaud]
- Make boundary a property for lazy computation. [Martin Raspaud]
- Create uptime on the fly if not provided. [Martin Raspaud]
- Add avoid_list feature. [Martin Raspaud]
- Fix ftp retrieval of aqua downlink schedule. [Martin Raspaud]

  Ftplib raises an error_perm sometime. It is now catched and handled
  correctly.
- Take into account the 'start' conf parameter. [Martin Raspaud]

  The 'start' config parameters aims at skipping the first passes of the
  schedule in order to avoid changing the next scheduled pass. It was
  unfortunately not being used at all. This patch fixes the code to the right
  behaviour.
- More debug info. [Martin Raspaud]
- Don't put whitespaces in plot filenames. [Martin Raspaud]
- Bugfixes and cleanup. [Martin Raspaud]
- Bugfix the bugfix. [Martin Raspaud]
- Merge pull request #3 from mraspaud/revert-2-develop. [Martin Raspaud]

  Revert "Change instrument from avhrr to avhrr/3"
- Revert "Change instrument from avhrr to avhrr/3" [Martin Raspaud]
- Merge pull request #2 from pnuu/develop. [Martin Raspaud]

  Change instrument from avhrr to avhrr/3
- Change instrument from avhrr to avhrr/3. [Panu Lahtinen]
- Merge pull request #1 from pnuu/simplified_platforms. [Martin Raspaud]

  Removed platform name to TLE translation
- Removed platform name to TLE translation. [Panu Lahtinen]
- Fix the case when last vertex of intersection was last vertex of
  polygon. [Martin Raspaud]
- Add setup.cfg for easy rpm generation. [Martin Raspaud]
- More spherical tests. [Martin Raspaud]
- Append tests to the test suite. [Martin Raspaud]
- Add a few test to spherical geometry. [Martin Raspaud]
- Add lons and lats to boundary init arguments. [Martin Raspaud]
- A None intersection now returns an area of 0. [Martin Raspaud]
- Update unittests to reflect structure changes. [Martin Raspaud]
- Put an example cfg in the base directory. [Martin Raspaud]
- Reorganizing. [Martin Raspaud]
- Shorter, more effective filenames for plots. [Martin Raspaud]
- Bugfix default xml location. [Martin Raspaud]
- Bugfix report function. [Martin Raspaud]
- Add reference area in plots. [Martin Raspaud]
- Add xml declarations for report mode. [Martin Raspaud]
- Add xml report mode. [Martin Raspaud]
- Make the graph option an input directory. [Martin Raspaud]
- Add option to generate pass plots. [Martin Raspaud]


v0.2.2 (2014-06-02)
-------------------
- Bump up version number. [Martin Raspaud]
- Sort passes to avoid conflicts. [Martin Raspaud]
- Add export method to graph. [Martin Raspaud]
- Fix backward compatibility issue with numpy. [Martin Raspaud]
- Refactorize, putting passes stuff in separate module. [Martin Raspaud]


v0.2.1 (2014-05-27)
-------------------

Fix
~~~
- Bugfix: wrong sorting of passes leaded to conflicting schedules.
  [Martin Raspaud]

Other
~~~~~
- Bump up version number. [Martin Raspaud]
- Make compare callable (as compare_scheds) [Martin Raspaud]
- Add the confirmation option to the compare script. [Martin Raspaud]
- Cleaning up. [Martin Raspaud]
- Add pykdtree to travis dependencies. [Martin Raspaud]


v0.2.0 (2014-05-20)
-------------------
- Bump up version number. [Martin Raspaud]
- Add option to compare the most recent requests to a confirmation.
  [Martin Raspaud]
- Save xml data to temporary file first. [Martin Raspaud]
- Refine station list. [Martin Raspaud]
- Add request/confirmation comparison. [Martin Raspaud]
- Remove dependency to scipy, and cleanup. [Martin Raspaud]
- Start the schedule a little before to make sure we don't start in the
  middle of a conflict. [Martin Raspaud]
- Added the glob dependency. [Martin Raspaud]
- If ftp can't be reached for aqua dumps, use cached data. [Martin
  Raspaud]
- Fix ftp export of xml file. [Martin Raspaud]
- Fix xml file ftp push. [Martin Raspaud]
- Add mail option to send errors by mail. [Martin Raspaud]
- Smallest passes allowed are 4 minutes long. [Martin Raspaud]
- Fix spherical intersection search. [Martin Raspaud]
- Run on euron1. [Martin Raspaud]
- Fix bug on intersection, where start of arc was the intersection.
  [Martin Raspaud]
- Added Bochum station. [Martin Raspaud]
- Added possibility to upload xmlfile to ftp. [Martin Raspaud]
- Add downloading of aqua dump times. [Martin Raspaud]
- Fix xml generation call. [Martin Raspaud]
- Add a few options in the config file. [Martin Raspaud]
- Use xml instead of lxml in the main xml generation function. [Martin
  Raspaud]
- Bugfix in installation requirements. [Martin Raspaud]
- Remove mpop from dependencies. [Martin Raspaud]
- Adding docs. [Martin Raspaud]
- Add atlas installation on travis. [Martin Raspaud]
- Added missing dependencies. [Martin Raspaud]
- Fixing travis. [Martin Raspaud]
- Renamed a few things to avoid -_ problems. [Martin Raspaud]
- Initial commit. [Martin Raspaud]
- Initial commit. [Martin Raspaud]
