Example
=======

How to generate a local reception schedule
------------------------------------------

This Bash-script shows an example how to use the schedule script::

	#!/usr/bin/bash
	. $HOME/.bash_profile
	bin=$HOME/.local/bin/schedule
	conf=$PYTROLL_BASE/etc-schedule
	logp=$PYTROLL_BASE/log/scheduler.log
	logc=$PYTROLL_BASE/log/create.log
	odir=$PYTROLL_BASE/schedules

	cfgfile=$1
	shift

	# report-mode with plots
	mode="-r -p"

	# min time-diff btwn passes
	delay="90"

	# create output for scisys
	sci="--scisys"

	# dont include aqua-dumps
	aqua="--no-aqua-dump"

	# write gv-files for dot
	graph="-g"

	# exit if no new tle
	(( `ls $PYTROLL_BASE/tle/*.tle 2>/dev/null|wc -l` > 0 )) || exit 0

	# catch newest TLE-file, remove others
	tle=`ls -t $PYTROLL_BASE/tle/*.tle|head -1`
	mv $tle $tle.txt
	rm -f $PYTROLL_BASE/tle/*.tle

	# settings for TLE-check
	satlst="NOAA 15|NOAA 18|NOAA 19|METOP-A|METOP-B|AQUA|SUOMI|TERRA "
	satcnt=8
	to_addr=pytroll@schedule

	# check if TLE-file is complete
	if (( `grep -P "$satlst" $tle.txt|tee .tle_grep|wc -l` != $satcnt )); then
	    exit 0
	else
	    tle="-t $tle.txt"
	fi

	# start schedule script
	$bin -v -l $logp -c $conf/$cfgfile $tle -d $delay -o $odir $graph $mode $sci $aqua $@


How to create instrument swath outlines in shapefile format for all passes scheduled
------------------------------------------------------------------------------------

To create instrument swath outlines in shapefile format one can run the
`create_shapefiles_from_schedule` script, e.g. like in the below example::

  create_shapefiles_from_schedule -s "Suomi NPP" NOAA-21 -t /data/tles/tle-202409060015.txt -x /data/schedules/2024-09-06-00-45-08-acquisition-schedule-confirmation-nrk.xml -o /data/pass_outlines

Then in the directory `/data/pass_outlines` a set of (rather small size)
shapefiles can be found, like the below one example::

  `viirs_Suomi-NPP_202409122257_202409122309_outline.cpg`
  `viirs_Suomi-NPP_202409122257_202409122309_outline.dbf`
  `viirs_Suomi-NPP_202409122257_202409122309_outline.prj`
  `viirs_Suomi-NPP_202409122257_202409122309_outline.shp`
  `viirs_Suomi-NPP_202409122257_202409122309_outline.shx`
