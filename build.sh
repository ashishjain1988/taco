#!/bin/bash

if !(command -v pyinstaller >/dev/null;) then
	echo \'pyinstaller\' is not installed on this machine. Exiting.
	exit 1
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
START_DIR=$PWD
python setup.py build
pyinstaller $DIR/refcomp_pyinstaller.spec
pyinstaller --onefile -n taco $DIR/taco/taco_run.py
mv $DIR/dist/taco "$START_DIR""/taco_run"
mv $DIR/dist/taco_refcomp "$START_DIR""/taco_refcomp"
