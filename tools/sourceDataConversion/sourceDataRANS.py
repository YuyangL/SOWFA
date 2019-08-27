# Function for processing SOWFA source term data averaged for RANS.
#
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 17 Nov 2016
# Modified by Yuyang Luan


# Figure out how many time directories there are within a directory, and put
# them in numerical order.
def getOutputTimes(dir):
  # Import necessary modules
  import os
  data = os.listdir(dir)
  outputTimesI = []
  outputTimes = []
  nTimes = len(data)
  ii = 0
  for i in range(nTimes):
    # unicode is replace by str in Python 3
    #  if (unicode(data[i][0]).isnumeric()):
     if (str(data[i][0]).isnumeric()):
        outputTimesI.append(data[i])
        ii = ii + 1

  nTimes = len(outputTimesI)

  outputTimesIndex = 0
  outputTimesSort = []
  for i in range(nTimes):
     outputTimesSort.append([i,float(outputTimesI[i])])

  outputTimesSort = sorted(outputTimesSort, key=lambda index: index[1])

  for i in range(nTimes):
     outputTimes.append(outputTimesI[outputTimesSort[i][0]])

  return nTimes, outputTimes


# Assemble a complete source history.
def assembleSourceHistory(inputDir):
  # Import necessary modules
  import numpy as np
  # Get the number of time directories and their names.
  [nTimes, outputTimes] = getOutputTimes(inputDir)
  # Initialize the big arrays.
  timeMomentumX = []
  timeMomentumY = []
  timeMomentumZ = []
  timeTemperature = []
  sourceMomentumX = []
  sourceMomentumY = []
  sourceMomentumZ = []
  sourceTemperature = []


  # Loop through the time directories and get the source information.
  for n in range(nTimes):
     sourceName = ['SourceUXHistory','SourceUYHistory','SourceUZHistory','SourceTHistory']
     for m in range(4):
        inputFile = inputDir + '/' + outputTimes[n] + '/' + sourceName[m]

        if (m == 0):
           [heightMomentum, timeMomentumXI, sourceMomentumXI] = readSourceHistoryFile(inputFile)

           if (n == 0):
               timeMomentumX = timeMomentumXI
               sourceMomentumX = sourceMomentumXI
           else:
               startTime = timeMomentumXI[0]

               l = len(timeMomentumX)

               if (timeMomentumX[l-1] > startTime):
                   indEnd = (np.where(timeMomentumX >= startTime))[0][0] - 1

               else:
                   indEnd = l-1

               timeMomentumX = np.append(timeMomentumX[0:indEnd],timeMomentumXI)
               sourceMomentumX = np.append(sourceMomentumX[0:indEnd][:],sourceMomentumXI,axis=0)

               timeMomentumXI = []
               sourceMomentumXI = []

        elif (m == 1):
           [heightMomentum, timeMomentumYI, sourceMomentumYI] = readSourceHistoryFile(inputFile)

           if (n == 0):
               timeMomentumY = timeMomentumYI
               sourceMomentumY = sourceMomentumYI
           else:
               startTime = timeMomentumYI[0]

               l = len(timeMomentumY)

               if (timeMomentumY[l-1] > startTime):
                   indEnd = (np.where(timeMomentumY >= startTime))[0][0] - 1

               else:
                   indEnd = l-1

               timeMomentumY = np.append(timeMomentumY[0:indEnd],timeMomentumYI)
               sourceMomentumY = np.append(sourceMomentumY[0:indEnd][:],sourceMomentumYI,axis=0)

               timeMomentumYI = []
               sourceMomentumYI = []

        elif (m == 2):
           [heightMomentum, timeMomentumZI, sourceMomentumZI] = readSourceHistoryFile(inputFile)

           if (n == 0):
               timeMomentumZ = timeMomentumZI
               sourceMomentumZ = sourceMomentumZI
           else:
               startTime = timeMomentumZI[0]

               l = len(timeMomentumZ)

               if (timeMomentumZ[l-1] > startTime):
                   indEnd = (np.where(timeMomentumZ >= startTime))[0][0] - 1

               else:
                   indEnd = l-1

               timeMomentumZ = np.append(timeMomentumZ[0:indEnd],timeMomentumZI)
               sourceMomentumZ = np.append(sourceMomentumZ[0:indEnd][:],sourceMomentumZI,axis=0)

               timeMomentumZI = []
               sourceMomentumZI = []

        elif (m == 3):
           [heightTemperature, timeTemperatureI, sourceTemperatureI] = readSourceHistoryFile(inputFile)

           if (n == 0):
               timeTemperature = timeTemperatureI
               sourceTemperature = sourceTemperatureI
           else:
               startTime = timeTemperatureI[0]

               l = len(timeTemperature)

               if (timeTemperature[l-1] > startTime):
                   indEnd = (np.where(timeTemperature >= startTime))[0][0] - 1

               else:
                   indEnd = l-1

               timeTemperature = np.append(timeTemperature[0:indEnd],timeTemperatureI)
               sourceTemperature = np.append(sourceTemperature[0:indEnd][:],sourceTemperatureI,axis=0)

               timeTemperatureI = []
               sourceTemperatureI = []

  return heightMomentum,heightTemperature,timeMomentumX,timeMomentumY,timeMomentumZ,timeTemperature,sourceMomentumX,sourceMomentumY,sourceMomentumZ,sourceTemperature


# Read a single source history file.
def readSourceHistoryFile(inputFile):
  import numpy as np

  # Open the file.
  fid = open(inputFile,'r')

  # Read the first line
  data = fid.readline()
  if (data[0] == 'H'):
     print('Variable with Height')

     # Get the source height information.
     heights = data
     iStart = heights.find(')')
     iEnd = heights.find('\n')
     heights = heights[iStart+2:iEnd-1]
     heights = np.array([float(s) for s in heights.split(' ')])

     # Close the source file.
     fid.close()

     # Read the source data and organize it.
    #  data = np.loadtxt(inputFile, dtype='string', skiprows=2)
    # 'string' doesn't work in Python 3
     data = np.loadtxt(inputFile, dtype=np.str, skiprows=2)
     time = np.transpose(np.array(data[:,0],dtype='float'))
     source = np.array(data[:,2:],dtype='float')

  elif (data[0] == 'T'):
     print('Constant with Height')

     heights = np.zeros(1);

     # Close the source file.
     fid.close()

     # Read the source data and organize it.
    #  data = np.loadtxt(inputFile, dtype='string', skiprows=1)
     data = np.loadtxt(inputFile,dtype=np.str,skiprows=1)
     time = np.transpose(np.array(data[:,0],dtype='float'))
     source = np.array(data[:,2:],dtype='float')

  return heights, time, source


def selectTimes(timeMomentumX, timeMomentumY, timeMomentumZ, timeTemperature, startTimes = (None, None, None, None), stopTimes = (None, None, None, None)):
    import numpy as np
    # Gather sources by momentum in x, y, z, and temperature
    timesAll = (timeMomentumX, timeMomentumY, timeMomentumZ, timeTemperature)
    # List treatment
    startTimes, stopTimes = list(startTimes), list(stopTimes)
    # Go through momentum x, y, z, and temperature
    iStart, iStop, startTimesReal, stopTimesReal = np.empty(4, dtype=np.int), np.empty(4, dtype=np.int), np.empty(4), np.empty(4)
    timesSelected = []
    for i in range(4):
        if startTimes[i] is None: startTimes[i] = timesAll[i][0]
        if stopTimes[i] is None: stopTimes[i] = timesAll[i][-1]
        # Bisection left to find actual starting and ending time and their indices
        (iStart[i], iStop[i]) = np.searchsorted(timesAll[i], (startTimes[i], stopTimes[i]))
        # If stopTime larger than any time, iStop = len(timesAll[i])
        iStop[i] = min(iStop[i], len(timesAll[i]) - 1)
        # Actual start and stop times
        startTimesReal[i], stopTimesReal[i] = timesAll[i][iStart[i]], timesAll[i][iStop[i]]
        # Append selected time list to a list
        timesSelected.append(timesAll[i][iStart[i]:iStop[i]])

    print('\nTime and index information extracted for ' + str(startTimesReal) + ' s - ' + str(stopTimesReal) + ' s for each momentum in x, y, z, and temperature')
    return timesSelected, startTimesReal, stopTimesReal, iStart, iStop


def calculateAverage(heightMomentum, heightTemperature,
 timeMomentumX, timeMomentumY, timeMomentumZ, timeTemperature,
 sourceMomentumX, sourceMomentumY, sourceMomentumZ, sourceTemperature, startTimes = (None, None, None, None), stopTimes = (None, None, None, None)):
    import numpy as np
    # Get selected times with given startTimes and stopTimes
    timesSelected, startTimesReal, stopTimesReal, iStart, iStop = selectTimes(timeMomentumX, timeMomentumY, timeMomentumZ, timeTemperature, startTimes, stopTimes)
    # Assign each selected time list to corresponding sources
    timeMomentumX_selected, timeMomentumY_selected, timeMomentumZ_selected, timeTemperature_selected = timesSelected[0], timesSelected[1], timesSelected[2], timesSelected[3]
    # Get 4 selected sources from selected times
    sourceMomentumX_selected, sourceMomentumY_selected, sourceMomentumZ_selected = sourceMomentumX[iStart[0]:iStop[0]],sourceMomentumY[iStart[1]:iStop[1]], sourceMomentumZ[iStart[2]:iStop[2]]
    sourceTemperature_selected = sourceTemperature[iStart[3]:iStop[3]]
    # Create a list of source*time according to heights
    sourceTimeMomentumX, sourceTimeMomentumY, sourceTimeMomentumZ = np.zeros_like(np.array(heightMomentum)), np.zeros_like(np.array(heightMomentum)), np.zeros_like(np.array(heightMomentum))
    sourceTimeTemperature = np.zeros_like(np.array(heightTemperature))
    # Go through each source in the order of momentum in x, y, z and calculate sum(srouce*time)
    for m in range(3):
        # For each source, go through selected time
        for n in range(len(timesSelected[m])):
            # For each time, go through each height and calculate sum(srouce*time)
            for i in range(len(heightMomentum)):
                if m == 0:
                    sourceTimeMomentumX[i] += sourceMomentumX_selected[n][i]*timeMomentumX_selected[n]
                elif m == 1:
                    sourceTimeMomentumY[i] += sourceMomentumY_selected[n][i]*timeMomentumY_selected[n]
                elif m == 2:
                    sourceTimeMomentumZ[i] += sourceMomentumZ_selected[n][i]*timeMomentumZ_selected[n]

    # Lastly, do the same for source temperature
    for n in range(len(timeTemperature_selected)):
        for i in range(len(heightTemperature)):
            sourceTimeTemperature[i] += sourceTemperature_selected[n][i]*timeTemperature_selected[n]

    # After going through each time, go through each height and do avg = sum(srouce*time)/sum(time)
    sourceMomentumX_mean, sourceMomentumY_mean, sourceMomentumZ_mean = np.empty_like(heightMomentum), np.empty_like(heightMomentum), np.empty_like(heightMomentum)
    sourceTemperature_mean = np.empty_like(heightTemperature)
    for i in range(len(heightMomentum)):
        sourceMomentumX_mean, sourceMomentumY_mean, sourceMomentumZ_mean = sourceTimeMomentumX/np.sum(timeMomentumX_selected), sourceTimeMomentumY/np.sum(timeMomentumY_selected), sourceTimeMomentumZ/np.sum(timeMomentumZ_selected)
        sourceTemperature_mean = sourceTimeTemperature/np.sum(timeTemperature_selected)

    return sourceMomentumX_mean, sourceMomentumY_mean, sourceMomentumZ_mean, sourceTemperature_mean, startTimesReal, stopTimesReal


# Write out the mean source file that will become SOWFA RANS input.
def writeMeanSourceForInput(fileName, heightMomentum,heightTemperature, sourceMomentumX_mean, sourceMomentumY_mean,sourceMomentumZ_mean, sourceTemperature_mean, startTimesReal, stopTimesReal):
    # Open the file.
    fid = open(fileName,'w')
    fid.write('// Sources averaged from {} s to {} s\n'.format(startTimesReal[0], stopTimesReal[0]))
    # Write the momentum source height list.
    fid.write('sourceHeightsMomentum\n')
    fid.write('(\n')
    for i in range(len(heightMomentum)):
        fid.write('    ' + str(heightMomentum[i]) + '\n')

    fid.write(');\n\n')
    # New placeholder time column for RANS
    timesNew = (0.0, 999999.9)
    # Write the x-momentum table
    fid.write('sourceTableMomentumX\n')
    fid.write('(\n')
    for n in range(len(timesNew)):
        textStr = '    (' + str(timesNew[n])
        for i in range(len(heightMomentum)):
            textStr = textStr + ' ' + str(sourceMomentumX_mean[i][0])

        textStr = textStr + ')\n'
        fid.write(textStr)

    fid.write(');\n')
    # Write the y-momentum table
    fid.write('sourceTableMomentumY\n')
    fid.write('(\n')
    for n in range(len(timesNew)):
        textStr = '    (' + str(timesNew[n])
        for i in range(len(heightMomentum)):
            textStr = textStr + ' ' + str(sourceMomentumY_mean[i][0])

        textStr = textStr + ')\n'
        fid.write(textStr)

    fid.write(');\n')
    # Write the z-momentum table
    fid.write('sourceTableMomentumZ\n')
    fid.write('(\n')
    for n in range(len(timesNew)):
        textStr = '    (' + str(timesNew[n])
        for i in range(len(heightMomentum)):
            textStr = textStr + ' ' + str(sourceMomentumZ_mean[i][0])

        textStr = textStr + ')\n'
        fid.write(textStr)

    fid.write(');\n')
    # Write the temperature source height list.
    fid.write('sourceHeightsTemperature\n')
    fid.write('(\n')
    for i in range(len(heightTemperature)):
        fid.write('    ' + str(heightTemperature[i]) + '\n')

    fid.write(');\n\n')
    # Write the temperature table
    fid.write('sourceTableTemperature\n')
    fid.write('(\n')
    for n in range(len(timesNew)):
        textStr = '    (' + str(timesNew[n])
        for i in range(len(heightTemperature)):
            textStr = textStr + ' ' + str(sourceTemperature_mean[n][i])

        textStr = textStr + ')\n'
        fid.write(textStr)

    fid.write(');\n')
    # Close the file.
    fid.close()
