import sourceDataRANS as sdRANS
import sys

# User input for the start and stop time of the temporal averaging of sources
print('sys = {}'.format(sys.argv))
startTime, stopTime = sys.argv[1], sys.argv[2]
# If start or stop time not provided, use None, i.e. no limit
startTime = None if startTime == '' else float(startTime)
stopTime = None if stopTime == '' else float(stopTime)

# Specify the directory where the source history files reside.
inputDir = './SourceHistory/'
outputFile = './sourcesRANS'

# Assemble the source information.
[heightMomentum,heightTemperature,
timeMomentumX,timeMomentumY,timeMomentumZ,timeTemperature,
sourceMomentumX,sourceMomentumY,sourceMomentumZ,sourceTemperature] = sdRANS.assembleSourceHistory(inputDir)

# Calculate temporal mean of sources
sourceMomentumX_mean, sourceMomentumY_mean, sourceMomentumZ_mean, sourceTemperature_mean = \
sdRANS.calculateAverage(heightMomentum, heightTemperature,
                        timeMomentumX, timeMomentumY, timeMomentumZ, timeTemperature,
                        sourceMomentumX, sourceMomentumY, sourceMomentumZ, sourceTemperature,  startTimes = (startTime,)*4, stopTimes = (stopTime,)*4)

# Write the mean source file for input to the windPlantADM.RANS solver
sdRANS.writeMeanSourceForInput(outputFile,
                            heightMomentum,
                            heightTemperature,
                            sourceMomentumX, sourceMomentumY,sourceMomentumZ, sourceTemperature)
