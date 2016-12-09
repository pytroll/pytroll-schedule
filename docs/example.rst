Example
=======

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
