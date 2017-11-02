Usage
=====

Usage of the schedule script::

	usage: schedule [-h] [-c CONFIG] [-t TLE] [-l LOG] [-m [MAIL [MAIL ...]]] [-v]
	                [--lat LAT] [--lon LON] [--alt ALT] [-f FORWARD]
	                [-s START_TIME] [-d DELAY] [-a AVOID] [--no-aqua-terra-dump]
	                [--multiproc] [-o OUTPUT_DIR] [-u OUTPUT_URL] [-x] [-r]
	                [--scisys] [-p] [-g]
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -c CONFIG, --config CONFIG
	                        configuration file to use
	  -t TLE, --tle TLE     tle file to use
	  -l LOG, --log LOG     File to log to (defaults to stdout)
	  -m [MAIL [MAIL ...]], --mail [MAIL [MAIL ...]]
	                        mail address(es) to send error messages to.
	  -v, --verbose         print debug messages too
	
	start-parameter:
	  (or set values in the configuration file)
	
	  --lat LAT             Latitude, degrees north
	  --lon LON             Longitude, degrees east
	  --alt ALT             Altitude, km
	  -f FORWARD, --forward FORWARD
	                        time ahead to compute the schedule
	  -s START_TIME, --start-time START_TIME
	                        start time of the schedule to compute
	  -d DELAY, --delay DELAY
	                        delay (in seconds) needed between two consecutive
	                        passes (60 seconds by default)
	
	special:
	  (additional parameter changing behaviour)
	
	  -a AVOID, --avoid AVOID
	                        xml request file with passes to avoid
	  --no-aqua-terra-dump  do not consider Aqua/Terra-dumps
	  --multiproc           use multiple parallel processes
	
	output:
	  (file pattern are taken from configuration file)
	
	  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
	                        where to put generated files
	  -u OUTPUT_URL, --output-url OUTPUT_URL
	                        URL where to put generated schedule file(s), otherwise
	                        use output-dir
	  -x, --xml             generate an xml request file (schedule)
	  -r, --report          generate an xml report file (schedule)
	  --scisys              generate a SCISYS schedule file
	  -p, --plot            generate plot images
	  -g, --graph           save graph info
