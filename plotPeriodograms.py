# plots periodograms

import matplotlib.pyplot as plt
import numpy as np
import os
import copy

def checkForTrueStatus(grams):
    count = 0
    for ID in grams.keys():
        stops = grams[ID]
        for stop in stops.keys():
            channels = stops[stop]
            for channel in channels.keys():
                if channels[channel]['status']:
                    count += 1

    print "true count: " + str(count)
    

def fillPeriodogramsDictStatus(path, filenames, notesfile):

    # creates dict with the format
    # dict['s1']['tp9']['status'] = True
    # dict['s1']['tp9']['data'] = np array periodogram
    # or
    # dict['s1']['tp9']['status'] = False
    # dict['s1']['tp9'] ['data'] = bad quality periodogram or []

    # might as well hardcode the nested dicts

    periodogramDict = {}
    
    # deep copy necessary
    # actually not if you bring into for scope below and keep creating the dict
    # actually that's what you should do Adam
    # actually maybe not, does this literally do the same thing but only create one?
    # is it faster? do we care? test dis out Adam my man
    keywords = ['all','bad','empty','missing']
    data = {'status' : False, 'data' : []}
    channels = {'tp9' : copy.deepcopy(data), 'tp10': copy.deepcopy(data),
                'af7': copy.deepcopy(data), 'af8' : copy.deepcopy(data)}
    stops = {'base' : copy.deepcopy(channels), 's1' : copy.deepcopy(channels),
             's2' : copy.deepcopy(channels), 's3' : copy.deepcopy(channels),
             's4' : copy.deepcopy(channels), 's5' : copy.deepcopy(channels),
             's6' : copy.deepcopy(channels)}


    
    # fill up dict with all data, updates status to True
    for filename in filenames:
        ID = filename.split('-')[0] # in case there is whitespace at end of line
        stop = filename.split('-')[1][:-4] # strips off .txt
        tp9,tp10,af7,af8 = getPeriodogramFromFile(path+filename)
        print len(periodogramDict.keys())
        print ID, stop

        if ID not in periodogramDict.keys():
            periodogramDict[ID] = copy.deepcopy(stops)

        periodogramDict[ID][stop]['tp9']['status'] = True
        periodogramDict[ID][stop]['tp9']['data'] = tp9
        periodogramDict[ID][stop]['tp10']['status'] = True
        periodogramDict[ID][stop]['tp10']['data'] = tp10
        periodogramDict[ID][stop]['af7']['status'] = True
        periodogramDict[ID][stop]['af7']['data'] = af7
        periodogramDict[ID][stop]['af8']['status'] = True
        periodogramDict[ID][stop]['af8']['data'] = af8
        
    checkForTrueStatus(periodogramDict)
    # now get status updates from notes
    with open(notesfile, 'r') as notes:
        lines = notes.readlines()

        print "number of files with notes: " + str(len(lines)/2.)

        for i in range(0,len(lines),2):

            ID_STOP = lines[i].strip() # in case there is whitespace at end of line
            ID, stop = ID_STOP.split('-')
            data = lines[i+1].strip()

            # we do not care about ID's in notes file if they are not in
            # periodogram files
            if ID in periodogramDict.keys():

                # keywords mean all channels are bad
                hit = False
                for key in keywords:
                    if key in data:
                        hit = True
                if hit:
                    for chan in channels.keys():
                        periodogramDict[ID][stop][chan]['status'] = False

                # now get specific channels
                for chan in channels.keys():
                    if chan in data:
                        periodogramDict[ID][stop][chan]['status'] = False
    checkForTrueStatus(periodogramDict)
    return periodogramDict
            
        

def getPeriodogramFromFile(filename):

    array = [[],[],[],[]]
    with open(filename,'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            line = line.split('\t')

            for i in range(4):
                array[i].append(float(line[i]))

    return array

# plots a single file on 4 subplots
# one participant at one stop, 4 electrodes
def plotFile(filename):

    pgrams = getPeriodogramFromFile(filename)

    freqs = np.linspace(55./len(pgrams[0]),55,len(pgrams[0]))

    for i in range(4):
        plt.subplot(2,2,i)
        plt.plot(freqs,pgrams[i])

    plt.show()

#plotFile('D:\\ORGANIZED\\PERIODOGRAMS\\999-s1.txt')


def plotChannel(channel, periodograms, freqs):

    stoplist = ['base','s1','s2','s3','s4','s5','s6']
    for ID in periodograms.keys():
        stops = periodograms[ID]
        for stop in stops.keys():
            channels = stops[stop]
            if channels[channel]['status']:
                plt.subplot(3,4,stoplist.index(stop)+1)
                plt.title(stop + channel)
                #plt.plot(freqs,20*np.log10(channels[channel]['data']))
                plt.plot(freqs, channels[channel]['data'])
                #print ID
                #print stop
                #print channels[channel]['data'][:5]

    plt.show()

def plotAverageChannel(channel, periodograms, freqs):

    stoplist = ['base','s1','s2','s3','s4','s5','s6']

    full_stops = [np.zeros(100000),np.zeros(100000),
                  np.zeros(100000),np.zeros(100000),
                  np.zeros(100000),np.zeros(100000),
                  np.zeros(100000)]
    counts = np.zeros(7)
    for ID in periodograms.keys():
        stops = periodograms[ID]
        for stop in stops.keys():
            channels = stops[stop]
            if channels[channel]['status']:
                print ID
                print stop
                full_stops[stoplist.index(stop)] += channels[channel]['data']
                counts[stoplist.index(stop)] += 1

    for x in range(7):
        plt.subplot(3,4,x+1)
        plt.title(stoplist[x] + channel)
        plt.plot(freqs, full_stops[x]/counts[x])
        #plt.plot(freqs,20*np.log10(full_stops[x]/counts[x]))

    plt.show()
    
                
def plotOneByOne(periodograms,freqs):
    for ID in periodograms.keys():
        stops = periodograms[ID]
        for stop in stops.keys():
            channels = stops[stop]
            count = 1
            for channel in channels.keys():
                if channels[channel]['status']:
                    plt.subplot(2,2,count)
                    plt.title(ID + stop + channel)
                    plt.plot(freqs,channels[channel]['data'])
                count+=1
            plt.show()

     
    
                
# main
#path = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\PERIODOGRAMS_BACKWARD\\"
#notes = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\NOTES_BACKWARD.txt"
#path = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\PERIODOGRAMS_TEST_FFTIGNORANT\\"
#notes = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\NOTES_TEST.txt"
path = "C:\\Users\\alzfranc\\Desktop\\ORGANIZED\\PERIODOGRAMS_FFTIG_ORIG_BACK\\"
#path = "C:\\Users\\alzfranc\\Desktop\\ORGANIZED\\PERIODOGRAMS_FFTIG_BACK\\"
notes = "C:\\Users\\alzfranc\\Desktop\\v2\\NOTES_BACKWARD.txt"
files = os.listdir(path)[0:150]
#b = [x[:x.index('-')] for x in a]
freqs = np.linspace(55./100000,55,100000)
fft_freqs = np.linspace(0, 110, 20001)

grams  = fillPeriodogramsDictStatus(path, files, notes)
#print grams['18168']['s6']['tp9']['status']

plotChannel('af7', grams, fft_freqs) #af7: 56

#plotAverageChannel('af7', grams, fft_freqs)
#plotOneByOne(grams,freqs)
#plotFile("C:\\Users\\alzfranc\\Desktop\\ORGANIZED\\PERIODOGRAMS_FFTIG_ORIG_BACK\\20703-s6.txt")
