#!/bin/bash

# 2023 © Vasileios Angelidakis <vasileios.angelidakis@fau.de>
# 2023 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>
# 2023 © Robert Caulk <rob.caulk@gmail.com>

# This script will run Tests 1-4 from Chung and Ooi (2011). It is standalone and should run on every linux system with yade installed.

# OUTPUTS:
# *.png and *.txt files per-job will be in ./outputData

# WARNING
# If for some reason this script execution is stopped with 'ctrl+z' before the end, then some jobs will stay behind and still use ressources, 
# don't forget 'pkill -9 yadedaily' after that (or whatever yade version it is) to really kill them

#export YADE="/path/to/my/own/yadeVersion"
#export YADE='yade'
export YADE='./yade-2021.01a'

export OMP_THREADS=1 # OpenMP threads, should be less than number of cores. Just 1 will suffice for these simple tests
export OMP_PROC_BIND=true # pin OMP threads to physical cores

# =========
#   TESTS  
# =========

# ----- TEST 1 -----
$YADE -j $OMP_THREADS -x Test1.py glass
$YADE -j $OMP_THREADS -x Test1.py limestone

# ----- TEST 2 -----
$YADE -j $OMP_THREADS -x Test2.py aluminium_alloy
$YADE -j $OMP_THREADS -x Test2.py magnesium_alloy

# ----- TEST 3 -----
for en in 1.0 0.8 0.6 0.4 0.2
do
	$YADE -j $OMP_THREADS -x Test3.py aluminium_oxide $en False # {True: en_es, False:betan_betas}
	$YADE -j $OMP_THREADS -x Test3.py cast_iron $en False
done

# ----- TEST 4 -----
for angle in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85
do
	$YADE -j $OMP_THREADS -x Test4.py aluminium_oxide $angle False # {True: en_es, False:betan_betas}
	$YADE -j $OMP_THREADS -x Test4.py aluminium_alloy $angle False
done


