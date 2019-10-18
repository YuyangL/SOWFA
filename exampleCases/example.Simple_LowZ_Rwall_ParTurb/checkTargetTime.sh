#!/usr/bin/env bash

# Requires 3 arguments: runNumber, targetEndTime and solver

# Check every $checkFreq s
checkFreq=60

# Sleep for 2 min when the solver initializes itself
sleep 120

# Itereative check until we forcoble kill the loop using break
while true
do
   # First read file from end to front, then look for "Time = ${2}"
   # -m 1 stop reading after 1st match
   # \b set a boundary so that ClockTime won't match
   if tac ./log.${1}.${3} | grep -m 1 "\bTime = ${2}"
   then
      # Change stopAt to writeNow in system/controlDict
      # \s means any number of tabs or spaces
      # .* means find pattern whatever is before .* and don't care what's after
      sed -i 's/stopAt\s.*/stopAt    writeNow;/' ./system/controlDict

      touch reachedTargetTime.stop

      # Break the loop
      break
   fi

   # Check every $checkFreq s
   sleep $checkFreq
done
