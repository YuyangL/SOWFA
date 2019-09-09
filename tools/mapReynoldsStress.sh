#!/usr/bin/env bash

# Requires 1 argument: [source directory]

module load openfoam/2.4.0

# Check for latestTime in the source directory and create the same in target directory
latestTime=`exec ls ${1} | sed 's/\([0-9]*\.[0-9]*\)/\1/g' | sort -n | tail -1`
mkdir $latestTime

# Make sure startTime at target's controlDict is ${latestTime} from source
sed -i "s/startTime\s.*/startTime    ${latestTime};/" ./system/controlDict

# Map from the latestTime from source
mapFields ${1} -fields '(uuPrime2)' -consistent -mapMethod 'mapNearest'

# Change startTime in target's controlDict back to 0
sed -i "s/startTime\s.*/startTime    0;/" ./system/controlDict
