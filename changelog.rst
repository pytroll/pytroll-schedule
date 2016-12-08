Changelog
=========

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


