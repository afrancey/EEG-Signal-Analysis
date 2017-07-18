# processes all the EEG files
# re-orders backward

#https://stackoverflow.com/questions/24702807/lomb-scargle-vs-fft-power-spectrum-crashes-with-evenly-spaced-data/
#lol

import os
from EEGProcessor import EEGSet
import numpy as np

nout = 4*100000 # was way too small at 1000!!
Fs = 220 # sampling rate
freqs = np.linspace(float(Fs)/nout, Fs, nout)
freqs = freqs[0:int(len(freqs)/4.)] # only looking at frequencies 0Hz to 55Hz
ang_freqs = 2*np.pi*np.linspace(float(Fs)/nout, Fs, nout)
ang_freqs = ang_freqs[0:int(len(ang_freqs)/4.)]

#pathToOrganized = 'F:\\'
pathToOrganized = "C:\\Users\\alzfranc\\Desktop\\"
pathToOrganized = "D:\\"

outputpath = pathToOrganized + 'ORGANIZED\\PERIODOGRAMS\\'

# forward
#forig = 'C:\\Users\\Adam Francey\\Desktop\\ORGANIZED\\data\\FORWARD\\EEG\\exported - filtered\\'
#freject = 'C:\\Users\\Adam Francey\\Desktop\\ORGANIZED\\data\\FORWARD\\EEG\\exported - post rejections\\'
#fnotes = 'notes'

forig = pathToOrganized + 'ORGANIZED\\data\\FORWARD\\EEG\\exported - filtered\\'
freject = pathToOrganized + 'ORGANIZED\\data\\FORWARD\\EEG\\exported - post rejections\\'
fnotes = 'notes'

# backward
#borig = 'C:\\Users\\Adam Francey\\Desktop\\ORGANIZED\\data\\BACKWARD\\EEG\\exported - filtered\\'
#breject = 'C:\\Users\\Adam Francey\\Desktop\\ORGANIZED\\data\\BACKWARD\\EEG\\exported - post rejections\\'
#bnotes = 'notes'

borig = pathToOrganized + 'ORGANIZED\\data\\BACKWARD\\EEG\\exported - filtered\\'
breject = pathToOrganized + 'ORGANIZED\\data\\BACKWARD\\EEG\\exported - post rejections\\'
bnotes = 'notes'

# new york
#norig = 'C:\\Users\\Adam Francey\\Desktop\\ORGANIZED\\data\\NEWYORK\\EEG\\exported - filtered\\'
#nreject = 'C:\\Users\\Adam Francey\\Desktop\\ORGANIZED\\data\\NEWYORK\\EEG\\exported - post rejections\\'
#nnotes = 'notes'

norig = pathToOrganized + 'ORGANIZED\\data\\NEWYORK\\EEG\\exported - filtered\\'
nreject = pathToOrganized + 'ORGANIZED\\data\\NEWYORK\\EEG\\exported - post rejections\\'
nnotes = 'notes'

forig = pathToOrganized + 'ORGANIZED\\data\\TEST\\exported - filtered\\'
freject = pathToOrganized + 'ORGANIZED\\data\\TEST\\exported - post rejections\\'

# dictionary of all filepaths
filepaths = {'forward':{'orig': forig, 'reject' : freject, 'notes':fnotes},
             'backward':{'orig': borig, 'reject' : breject, 'notes':bnotes},
             'newyork': {'orig': norig, 'reject' : nreject, 'notes':nnotes}
            }


def getListOfIDs(filenames):

    IDs = []
    for name in filenames:
        IDs.append(name[:name.rindex('-')]) # everything before last (second) dash
        
    return set(IDs)


# forward to test

stops = ['base','s1','s2','s3','s4','s5','s6']

paths = filepaths['newyork']

origfiles = os.listdir(paths['orig'])
rejectfiles = os.listdir(paths['reject'])

IDs = getListOfIDs(origfiles)

setCount = 0.
bugCount = 0.
missingCount = 0.
valueCount = 0.
zeroDivision = 0.
errorCount=0.


stringToWrite = ""

import time
startTime = time.time()
for ID in IDs:
    #print ID

    try:
        UUID = ID.split('-')[1]
    except IndexError:
        print "MINOR ERROR: NO TIMESTAMP IN ID"
        UUID = ID

    for stop in stops:
        #print stop

        orig = paths['orig'] + ID + '-' + stop + '_filtered.txt'
        reject = paths['reject'] + UUID + '-' + stop + '_rejected.txt'
        events = paths['reject'] + UUID + '-' + stop + '_eventSet.txt'

        print orig
        eset =  EEGSet(orig,reject,events)
        setCount += 1

        if eset.okayToProcess:
            try:
                pgrams, bandpowers, relative = eset.process(freqs)
                with open(outputpath + UUID + "-" + stop + ".txt",'w') as f:

                    for sample in range(len(pgrams[0])):
                        f.write(str(pgrams[0][sample]) + '\t' + str(pgrams[1][sample]) + '\t'
                            + str(pgrams[2][sample]) + '\t' + str(pgrams[3][sample]) + '\n')
                                
                                
            except ZeroDivisionError:
                print "ZERO DIVISION ERROR"
                zeroDivision+=1
        else:
            print "ERROR: " + eset.error
            errorCount+=1



        # legacy after setCount +=1 --------
            
##        if eset.indicatorArray == False:
##            bugCount += 1
##            print "PERCENT AFFECTED: " + str(100*bugCount/setCount)
##            stringToWrite += orig + "\n"
##        else:
##            if eset.TEST_indicator(eset.originalSet,eset.rejectSet, eset.indicatorArray) and eset.rejectSet != "file does not exist":
##                try:
##                    pgrams, bandpowers, relative = eset.process(freqs)
##                    with open(outputpath + UUID + "-" + stop + ".txt",'w') as f:
##
##                        for sample in range(len(pgrams[0])):
##                            f.write(str(pgrams[0][sample]) + '\t' + str(pgrams[1][sample]) + '\t'
##                                    + str(pgrams[2][sample]) + '\t' + str(pgrams[3][sample]) + '\n')
##                                
##                                
##                except ZeroDivisionError:
##                    print "ZERO DIVISION ERROR"
##                    zeroDivision+=1

        # legacy-------
                    

    
    

print IDs
print len(IDs)

print stringToWrite
print setCount
print bugCount
print missingCount
print valueCount
print zeroDivision
print errorCount
print "time: " + str(time.time() - startTime)
print freqs

with open("failures_NAN.txt",'w') as f:
    f.write(stringToWrite)
    


# gets a list of ids in a set
# expected filenames named such as 1442676390786-17635-s3_filtered.txt'

        

    









# to do... make this read csv and also make the csv file
def get_channelsToReject_CSV(direction):

    if direction == 'forward':
        notesfilename = 'C:/Users/Adam Francey/Desktop/eeglab/to export/EEG Forward new data/NOTES_Forward.txt'
    elif direction == 'backward':
        notesfilename = 'C:/Users/Adam Francey/Desktop/eeglab/to export/EEG Reversed New Data/NOTES_Reversed.txt'
    elif direction =='new york':
        notesfilename = 'C:/Users/Adam Francey/Desktop/eeglab/to export/new york new data/NOTES_NEWYORK.txt'

    suspects_tp9 = []
    suspects_tp10 = []
    suspects_af7 = []
    suspects_af8 = []

    with open(notesfilename,'r') as notes:

        lines = notes.readlines()

        for line in lines:
            line = line.strip('\n').split('\t')

            id_stop = line[0]
            note = line[1]

            regularNote = False

            if 'tp9' in note:
                suspects_tp9.append(id_stop) 
                regularNote = True

            if 'tp10' in note:
                suspects_tp10.append(id_stop)
                regularNote = True

            if 'af7' in note:
                suspects_af7.append(id_stop) 
                regularNote = True

            if 'af8' in note:
                suspects_af8.append(id_stop)
                regularNote = True

            if 'empty' in note or 'missing' in note or 'all' in note or 'unable' in note:

                # '-unable to read "not a binary MAT-file"''
                # no channels named - all to be rejected
                suspects_tp9.append(id_stop) 
                suspects_tp10.append(id_stop)
                suspects_af7.append(id_stop) 
                suspects_af8.append(id_stop)
                regularNote = True

            if regularNote == False:
                print note




    return suspects_tp9, suspects_tp10, suspects_af7,suspects_af8
