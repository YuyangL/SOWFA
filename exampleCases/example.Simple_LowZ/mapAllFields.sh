#!/usr/bin/env bash

# Requires 1 argument: [source directory] [time in source]

# module load openfoam/2.4.0

targetDir=$PWD

# Check for specified time in the source directory and create the same in target directory
if [ $2 == 'latestTime' ]
    then
    latestTime=`exec ls ${1} | sed 's/\([0-9]*\.[0-9]*\)/\1/g' | sort -n | tail -1`
    mkdir $latestTime
else
    $latestTime=$2
    mkdir $2
fi

# Make sure startTime at target's controlDict is the specified time $2 from source
sed -i "s/startTime\s.*/startTime    ${2};/" ./system/controlDict
# Map from the specified time $2 from source in $1 dir
mapFields $1 -sourceTime $2  -fields '(U T p_rgh qwall Rwall k epsilonSGS nuSgs kappat)' -consistent -mapMethod 'mapNearest'

# Change startTime in target's controlDict back to 0
sed -i "s/startTime\s.*/startTime    0;/" ./system/controlDict

# Copy al fields from ${latestTime} dir in RANS case to 0/ of RANS case and remove $latestTime folder
rsync -avh $latestTime/* 0/
rm -rf $latestTime

# Lastly, rename k*, epsilon*, and nu* names to k, epsilon and nut
mv 0/k* 0/k.gz
mv 0/epsilon* 0/epsilon.gz
mv 0/nu* 0/nut.gz
