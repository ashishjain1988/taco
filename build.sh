#!/bin/bash
set -e
set -o pipefail

if !(command -v pyinstaller >/dev/null;) then
	echo \'pyinstaller\' is not installed on this machine. Exiting.
	exit 1
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
START_DIR=$PWD
python setup.py build
cp -r $DIR/share build/lib*/taco
cp -r $DIR/refcomp_pyinstaller.spec build/lib*/taco/
cd build/lib*/taco
pyinstaller refcomp_pyinstaller.spec
pyinstaller --onefile taco_run.py
mv dist/taco_run "$START_DIR""/taco_run"
mv dist/taco_refcomp "$START_DIR""/taco_refcomp"
