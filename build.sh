#!/bin/bash

START_DIR=$PWD
python setup.py build
cd build
cd lib*
cd taco
pyinstaller --onefile taco_run.py
cd dist
mv taco_run "$START_DIR""/taco_run"
