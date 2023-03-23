#! /bin/bash

orgdir=`pwd`
exertdir=`dirname $(realpath $0)`"/.."

cd $exertdir #EXERT

git clone --recursive https://github.com/brentp/python-giggle pyGiggle
cd pyGiggle
python setup.py test
python setup.py install