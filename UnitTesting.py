# Unit Testing for EEG Processor
import numpy as np
from EEGProcessor import EEGSet
import os

import pathlib

class Empatica():

    def __init__(self, folderpath, datalist):

        # make the folder

        self.folderpath = folderpath
        self.datalist = datalist

        try:
            pathlib.Path(folderpath).mkdir(parents=True, exist_ok=True)
        except:
            print("cannot create existing folder")

        # make other files
        otherfilenames = ['ACC.csv', 'BVP.csv', 'HR.csv', 'IBI.csv', 'info.txt',
                          'tags.csv', 'TEMP.csv']

        # literally create the files and do nothing else
        for filename in otherfilenames:
            with open(folderpath + "/" + filename, 'w') as f:
                pass

        # make the datafile

        with open(folderpath + "/" + 'EDA.csv', 'w') as f:
            f.write('1583285975\n4.0' + '\n'.join(datalist))

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

        print("Output string for this file:")
        print(eset.output_string)
        
    

                                     
#print("time: " + str(time.time() - startTime))
#print(freqs)

# TEST EDA
# make participant file

print("##################################################")
print(" MAKING EDA FILES")

try:
    os.mkdir("C:/Users/alzfr/Desktop/testEDA/UnitTesting files/raw files/const10/")
except:
    print("cannot create existing folder")
    
with open("C:/Users/alzfr/Desktop/testEDA/UnitTesting files/raw files/const10/EDA.csv", "w") as f:
    f.write("1524483828.23397\n4.0000")
    values = [str(10) for x in range(1000)]
    stringToWrite = ""
    for v in values:
        stringToWrite+="\n" + v
    f.write(stringToWrite)


boundaries = ["const10,0 100 200 250 333 448 882 932 "]
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

for filename in os.listdir(inputpath):

    print("Calculating (mean, slope) for file: " + filename)
    eset = EEGSet(inputpath + filename, boundaryfilepath, "EDA")
    mean, slope = eset.process()
    print("Should see mean = 10 and slope = 0")
    print(str(mean) + ", " + str(slope))

    print("Output string for this file: ")
    print(eset.output_string)


# create Empatica files to test whole pipeline
emp1data = [str(10) for x in range(180*4)] # 3 mins garbage
emp1data += [str(5) for x in range(20*4)] # 25 secs good
emp1data += [str(0) for x in range(4*4)] # 4 secs artifact
emp1data += [str(5) for x in range(30*4)] # 30 secs good
emp1data += [str(0) for x in range(1*4)] # 1 secs artifact
emp1data += [str(5) for x in range(15*4)] # 15 secs good
emp1 = Empatica('C:/Users/alzfr/Desktop/testEDA/UnitTesting files/testfiles/emp1', emp1data)

emp2data = [str(10) for x in range(180*4)] # 3 mins garbage
emp2data += [str(0) for x in range(20*4)] # 25 secs art
emp2data += [str(5) for x in range(4*4)] # 4 secs good
emp2data += [str(0) for x in range(30*4)] # 30 secs art
emp2data += [str(5) for x in range(1*4)] # 1 secs good
emp2data += [str(0) for x in range(15*4)] # 15 secs art
emp2 = Empatica('C:/Users/alzfr/Desktop/testEDA/UnitTesting files/testfiles/emp2', emp2data)

# after inspecting files
inputpath = "C:/Users/alzfr/Desktop/testEDA/UnitTesting files/testfiles/" # folder which contains EEG files
boundaryfilepath = "C:/Users/alzfr/Desktop/testEDA/UnitTesting files/bounds_1583451288194.csv" # path to boundaries
outputfilepath = "C:/Users/alzfr/Desktop/testEDA/UnitTesting files/inspectiontestoutput.csv"

stringtowrite = ''
for filename in os.listdir(inputpath):

    if 'config' not in filename:
        print("Calculating (mean, slope) for file: " + filename)
        eset = EEGSet(inputpath + filename, boundaryfilepath, "EDA")
        mean, slope = eset.process()
        print("Should see mean = 5 and slope = 0")
        print(str(mean) + ", " + str(slope))

        print("Output string for this file: ")
        print(eset.output_string)
        stringtowrite += eset.output_string + "\n"


with open(outputfilepath, 'w') as f:
    f.write(stringtowrite)


    



