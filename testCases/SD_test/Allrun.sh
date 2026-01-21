#!/bin/bash
#
#-------------------------------------------------#
# Contributor: Mario Javier Rincon                #
# Updated on:  10 December 2020                   #
#-------------------------------------------------#
# Topic:       Optimisation Flow Meters           #
# OpenFOAM:    v2006                              #
#-------------------------------------------------#
# Institution: Aarhus University / Kamstrup       #
# Email:       mjrp@eng.au.dk                     #
#-------------------------------------------------#
#
#------------------------------------------------------------------------------
# Source OpenFOAM run functions
# .$WM_PROJECT_DIR/bin/tools/RunFunctions

#------------------------------------------------------------------------------
#cd ${0%/*} || exit 1 # Run from this directory, otherwise exit
# clear

touch STARTFILE

rm -r logs

foamListTimes -rm

cp -r 0.orig 0

mkdir logs
#------------------------------------------------------------------------------
echo -e "   - Running simpleFoam"
simpleFoam > logs/simulation

#------------------------------------------------------------------------------
echo -e "   - Running post processing"
python post_process.py

foamLog logs/simulation

python plot_residuals.py

#------------------------------------------------------------------------------
echo -e "   - Simulation done :)"

touch DONEFILE
#------------------------------------------------------------------------------
