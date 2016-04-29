#!/bin/bash

if !(command -v pyinstaller >/dev/null;) then
	echo \'pyinstaller\' is not installed on this machine. Exiting.
	exit 1
fi

START_DIR=$PWD
python setup.py build
cd build
cd lib*
cd taco
pyinstaller --onefile -n taco_run taco_run.py
cd dist
mv taco_run "$START_DIR""/taco_run"
