# processes all the EEG files
# re-orders backward

#https://stackoverflow.com/questions/24702807/lomb-scargle-vs-fft-power-spectrum-crashes-with-evenly-spaced-data/
#lol

import os
from EEGProcessor import EEGSet
import numpy as np

nout = 4*100000 # What is proper size? Any point in oversampling?
Fs = 220 # sampling rate
freqs = np.linspace(float(Fs)/nout, Fs, nout)
freqs = freqs[0:int(len(freqs)/4.)] # only looking at frequencies 0Hz to 55Hz
ang_freqs = 2*np.pi*freqs

inputpath = "" # folder which contains EEG files
boundaryfile = "" # path to boundaries
outputpath = ""

stringToWrite = ""

import time
startTime = time.time()

for filename in os.listdir(inputpath):
    if "EEG" in filename:
        eset = EEGSet(inputpath + "/" + filename)
        pgrams, bandpowers, relative = eset.process(freqs)
                    
    
print("time: " + str(time.time() - startTime))
print(freqs)
