# processes all the EEG files
# re-orders backward

#https://stackoverflow.com/questions/24702807/lomb-scargle-vs-fft-power-spectrum-crashes-with-evenly-spaced-data/
#lol

import os
from EEGProcessor import EEGSet
import numpy as np



nout = 4*100000

# nout = number of evenly-spaced frequencies to evaluate
# Choosing appropriate nout
# If
# N = num samples,
# df = spacing between frequencies in periodogram
# dt = time between samples
# then we require
# df >= 1/(N*dt)
# source: Jacob T. VanderPlas, "Understanding the Lomb-Scargle Periodogram", pg 14
#           https://arxiv.org/pdf/1703.09824.pdf
                
Fs = 220 # sampling rate
freqs = np.linspace(float(Fs)/nout, Fs, nout)
freqs = freqs[0:int(len(freqs)/4.)] # only looking at frequencies 0Hz to 55Hz
ang_freqs = 2*np.pi*freqs

inputpath = "" # folder which contains EEG files
boundaryfilepath = "" # path to boundaries
outputpath = ""

stringToWrite = ""

import time
startTime = time.time()

for filename in os.listdir(inputpath):
    if "EEG" in filename:

        print("Calculating periodograms for file: " + filename)
        eset = EEGSet(inputpath + "/" + filename, boundaryfilepath)
        pgrams, bandpowers, relative = eset.process(freqs)

        stringToWrite+= filename + ","

        #relative[i][j] = band power at channel i band j
        
        for band in range(0,5):
            stringToWrite+= ",".join([str(relative[channel][band]) for channel in range(0,4)] + ","
        
                    
    
print("time: " + str(time.time() - startTime))
print(freqs)
