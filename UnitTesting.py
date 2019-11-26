# Unit Testing for EEG Processor
import numpy as np
from EEGProcessor import EEGSet
import os

# make files

# file number 1:
# 20 secs 5 Hz, 20 secs 10 Hz, 20 secs 20 Hz
# first three channels have boundaries covering two sections each
# last channel has boundaries spread out

time_points = np.linspace(0, 60 - 1/220, 60*220)
data_points = np.append(np.sin(2*np.pi*5*time_points[0:20*220]), np.append(np.sin(2*np.pi*10*time_points[20*220:40*220]),np.sin(2*np.pi*20*time_points[40*220:60*220])))

# make participant file
with open("C:/Users/alzfr/Desktop/testLombscargle/UnitTesting files/filtered files/EEGsine,sine.txt", "w") as f:
    f.write("index\ttp9\ttp10\tfp1\tfp2\n")
    stringtowrite = ""
    for i in range(len(time_points)):
        stringtowrite+=str(i) + "\t" + str(data_points[i]) + "\t" + str(data_points[i]) + "\t" + str(data_points[i]) + "\t" + str(data_points[i]) + "\n"
    f.write(stringtowrite)
        

# file 2: clean set

boundaries = ["EEGsine,sine.txt,0 8800 ,0 4400 8800 13200 ,4400 13200 ,400 600 2746 3422 7882 8999 9032 10000 "]
# make boundaries file
with open("C:/Users/alzfr/Desktop/testLombscargle/UnitTesting files/boundstest.csv", "w") as f:
    for i in boundaries:
        f.write(i + "\n")


# TEST EEG
inputpath = "C:/Users/alzfr/Desktop/testLombscargle/UnitTesting files/filtered files/" # folder which contains EEG files
boundaryfilepath = "C:/Users/alzfr/Desktop/testLombscargle/UnitTesting files/boundstest.csv" # path to boundaries
outputfilepath = "C:/Users/alzfr/Desktop/testLombscargle/UnitTesting files/output.csv"

stringToWrite = ""

import time
startTime = time.time()

for filename in os.listdir(inputpath):
    if "EEG" in filename:

        print("Calculating periodograms for file: " + filename)
        eset = EEGSet(inputpath + filename, boundaryfilepath, "EEG")
        pgrams, bandpowers, relative = eset.process()

        # Test output
        print("Band Order: delta, theta, alpha, beta")
        print("EEGsine")
        print("Channel 1 Bandpowers: Expect to see large beta")
        print(bandpowers[0])
        print("Channel 2: Expect to see large alpha")
        print(bandpowers[1])
        print("Channel 3: Expect to see large theta")
        print(bandpowers[2])
        print("Channel 4: Expect to see a mix")
        print(bandpowers[3])
        
    

                                     
#print("time: " + str(time.time() - startTime))
#print(freqs)

# TEST EDA
# make participant file

print("##################################################")
print(" MAKING EDA FILES")
with open("C:/Users/alzfr/Desktop/testEDA/UnitTesting files/raw files/const10.txt", "w") as f:
    f.write("1524483828.23397\n4.0000")
    values = [str(10) for x in range(1000)]
    stringToWrite = ""
    for v in values:
        stringtowrite+="\n" + v
    f.write(stringtowrite)


boundaries = ["const10.txt,0 100 200 250 333 448 882 932 "]
# make boundaries file
with open("C:/Users/alzfr/Desktop/testEDA/UnitTesting files/boundstest.csv", "w") as f:
    for i in boundaries:
        f.write(i + "\n")

inputpath = "C:/Users/alzfr/Desktop/testEDA/UnitTesting files/raw files/" # folder which contains EEG files
boundaryfilepath = "C:/Users/alzfr/Desktop/testEDA/UnitTesting files/boundstest.csv" # path to boundaries
outputfilepath = "C:/Users/alzfr/Desktop/testEDA/UnitTesting files/output.csv"

stringToWrite = ""

import time
startTime = time.time()

print("...")
print(" RESULTS")
print("")

##for filename in os.listdir(inputpath):
##    if "EDA" in filename:
##
##        print("Calculating periodograms for file: " + filename)
##        eset = EEGSet(inputpath + filename, boundaryfilepath)
##        pgrams, bandpowers, relative = eset.process()
##
##        # Test output
##        print("Band Order: delta, theta, alpha, beta")
##        print("EEGsine")
##        print("Channel 1 Bandpowers: Expect to see large beta")
##        print(bandpowers[0])
##        print("Channel 2: Expect to see large alpha")
##        print(bandpowers[1])
##        print("Channel 3: Expect to see large theta")
##        print(bandpowers[2])
##        print("Channel 4: Expect to see a mix")
##        print(bandpowers[3])

