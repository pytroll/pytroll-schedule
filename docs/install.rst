Installation
============
	mkdir $HOME/PyTROLL
	cd !$
	mkdir data-in data-out etc

	#dnf install pyshp python-configobj numpy numpy-f2py scipy python-numexpr \
	#python-pillow proj proj-epsg proj-nad pyproj python2-matplotlib python-basemap
	#
	wget https://github.com/pytroll/mipp/archive/master.zip -O mipp-master.zip
	wget https://github.com/pytroll/mpop/archive/pre-master.zip -O mpop-pre-master.zip
	wget https://github.com/pytroll/pyresample/archive/master.zip -O pyresample-master.zip
	wget https://github.com/adybbroe/python-geotiepoints/archive/master.zip -O python-geotiepoints-master.zip
	wget https://github.com/pytroll/pycoast/archive/master.zip -O pycoast-master.zip
	wget https://github.com/pytroll/pyorbital/archive/master.zip -O pyorbital-master.zip
	wget https://github.com/pytroll/trollcast/archive/master.zip -O trollcast-master.zip
	wget https://github.com/pytroll/pytroll-schedule/archive/master.zip -O pytroll-schedule-master.zip
	wget https://github.com/pytroll/trollduction/archive/master.zip -O trollduction-master.zip


	for ff in *.zip ; do
	   unzip -u $ff
	   cd $(basename $ff .zip)
	   python setup.py build
	   python setup.py install --user
	done

	cat <<EE
	export PPP_CONFIG_DIR=$HOME/PyTROLL/etc
	export XRIT_DECOMPRESS_OUTDIR=$HOME/PyTROLL/data-out
	EE > $HOME/.pytroll.rc
