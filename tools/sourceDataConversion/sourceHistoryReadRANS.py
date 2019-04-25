import sourceDataRANS as sdRANS
import sys
from warnings import warn

# User input for the start and stop time of the temporal averaging of sources
print('sys = {}'.format(sys.argv))
try:
    startTime = sys.argv[1]
except:
    warn('\nNo start time provided as 1st argument. Assuming beginning time for averaging!\n', stacklevel = 2)
    # If start or stop time not provided, use None, i.e. no limit
    startTime = None

try:
    stopTime = sys.argv[2]
except IndexError:
    warn('\nNo stop time provided as 2nd argument. Assuming non-stop for averaging!\n', stacklevel = 2)
    stopTime = None

# Specify the directory where the source history files reside.
inputDir = './SourceHistory/'
outputFile = './sourcesRANS'

# Assemble the source information.
[heightMomentum,heightTemperature,
timeMomentumX,timeMomentumY,timeMomentumZ,timeTemperature,
sourceMomentumX,sourceMomentumY,sourceMomentumZ,sourceTemperature] = sdRANS.assembleSourceHistory(inputDir)

# Calculate temporal mean of sources
sourceMomentumX_mean, sourceMomentumY_mean, sourceMomentumZ_mean, sourceTemperature_mean, startTimesReal, stopTimesReal = \
sdRANS.calculateAverage(heightMomentum, heightTemperature,
                        timeMomentumX, timeMomentumY, timeMomentumZ, timeTemperature,
                        sourceMomentumX, sourceMomentumY, sourceMomentumZ, sourceTemperature,  startTimes = (startTime,)*4, stopTimes = (stopTime,)*4)

# Write the mean source file for input to the windPlantADM.RANS solver
sdRANS.writeMeanSourceForInput(outputFile,
                            heightMomentum,
                            heightTemperature,
                            sourceMomentumX, sourceMomentumY,sourceMomentumZ, sourceTemperature,
                            startTimesReal, stopTimesReal)

print('\nsourcesRANS saved after averaging from {} s to {} s. Note that the file needs to be renamed to "sources"'.format(startTimesReal[0], stopTimesReal[0]))
