# processes all the EEG files
# re-orders backward

#https://stackoverflow.com/questions/24702807/lomb-scargle-vs-fft-power-spectrum-crashes-with-evenly-spaced-data/
#lol

import os
from EEGProcessor import EEGSet
import numpy as np

inputpath = "C:/Users/alzfr/Desktop/testLombscargle/filtered (1,30) order 3 data/" # folder which contains EEG files
boundaryfilepath = "C:/Users/alzfr/Desktop/testLombscargle/inspected/combined.csv" # path to boundaries
outputfilepath = "C:/Users/alzfr/Desktop/testLombscargle/output.csv"

stringToWrite = ""

import time
startTime = time.time()

for filename in os.listdir(inputpath):
    if "EEG" in filename:

        print("Calculating periodograms for file: " + filename)
        eset = EEGSet(inputpath + filename, boundaryfilepath)
        pgrams, bandpowers, relative = eset.process(freqs)

        stringToWrite+= filename + ","

        #relative[i][j] = band power at channel i band j
        
        for band in range(0,5):
            print("help")
            
            stringToWrite+= ",".join([str(relative[channel][band]) for channel in range(0,4)]) + ","
                                     
#print("time: " + str(time.time() - startTime))
#print(freqs)
