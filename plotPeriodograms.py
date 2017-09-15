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

    freqs = np.linspace(55./len(pgrams[0]),55,len(pgrams[0])) # IF PLOTTING FFT PERIODOGRAM, FREQUENIES ACTUALLY GO TO 110!
    channels = ['tp9', 'tp10', 'af7', 'af8']
    for i in range(1,5):
        plt.subplot(2,2,i)
        plt.title(channels[i-1])
        plt.plot(freqs,pgrams[i-1])

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

def plotSignalFile(filename):

    array = [[],[],[],[]]
    with open(filename,'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            line = line.strip()
            line = line.split('\t')

            for i in range(4):
                array[i].append(float(line[i]))

    for i in range(1,5):
        plt.subplot(2,2,i)
        plt.plot(array[i-1][:1000])

    plt.show()

def plotChannelFromFile(filename, channel):

    # filename: filename including path
    # channel: channel
    # ptype: type of periodogram: lomb or fft

    pgrams = getPeriodogramFromFile(filename)
    freqs = np.linspace(55./len(pgrams[0]),55,len(pgrams[0])) # IF PLOTTING FFT PERIODOGRAM, FREQUENIES ACTUALLY GO TO 110!
        
    channels = ['tp9', 'tp10', 'af7', 'af8']
    plt.plot(freqs,pgrams[channels.index(channel)])
    #plt.plot(pgrams[channels.index(channel)])
    #plt.show()


def save_band_powers(freqs, listof_intervals, listof_band_names, periodograms, amp_spec_dens = True):

    # mutates periodograms, adds the band power list to channel dict
    # all periodograms were created as amplitude spectral densities
    # amp_spec_dens = True reverses this tranform to get power

    # first translate frequency intervals to indicies
    maxFreq = freqs[len(freqs)-1]
    index_intervals = []
    for band in listof_intervals:
        start = band[0]
        end = band[1]
        start_i = int((start/maxFreq)*len(freqs))
        end_i = int((end/maxFreq)*len(freqs))
        index_intervals.append([start_i,end_i])
        
    # calculate band powers
    for ID in periodograms.keys():
        stops = periodograms[ID]
        for stop in stops.keys():
            channels = stops[stop]
            for channel in channels.keys():
                if channels[channel]['status']:

                    pgram = np.array(channels[channel]['data'])
                    if amp_spec_dens:
                        pgram = (pgram**2)*len(pgram)/4

                    total_power = np.sum(pgram)

                    # get band power by summing indicies
                    band_powers = []
                    for band in index_intervals:
                        start = band[0]
                        end = band[1]
                        band_power = np.sum(pgram[start:end])/total_power
                        band_powers.append(band_power)

                    # add it to periodogram dict (mutate)
                    channels[channel]['band powers'] = band_powers

    # save file
    stops_ordered = ['base', 's1','s2','s3','s4','s5','s6']
    channels_ordered = ['tp9','tp10','af7','af8']

    with open("BAND_POWERS.csv", 'w') as f:

        # write header
        header_to_write = 'ID'
        for s in stops_ordered:
            for c in channels_ordered:
                for b in listof_band_names:
                    header_to_write += ',' + s + c + b
        header_to_write += '\n'
        f.write(header_to_write)
                    
        
        # write line for every participant
        for ID in periodograms.keys():
            line_to_write = ID
            stops = periodograms[ID]
            for stop in stops_ordered:
                channels = stops[stop]
                for channel in channels_ordered:
                    if channels[channel]['status']:
                        band_powers = channels[channel]['band powers']
                        for band_power in band_powers:
                            line_to_write += ',' + str(band_power)
                    else:
                        # no band powers here
                        for i in range(len(listof_band_names)):
                            line_to_write += ',NAN'
            line_to_write += '\n'
            f.write(line_to_write)
            
    
                
# main
#path = 'D:\\ORGANIZED\\PERIODOGRAMS_TEST_3windows_halfoverlap\\'
#path = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\PERIODOGRAMS_BACKWARD\\"
#notes = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\NOTES_BACKWARD.txt"
#path = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\PERIODOGRAMS_TEST_FFTIGNORANT\\"
#path = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\PERIODOGRAMS_TEST\\"
#path = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\PERIODOGRAMS_TEST_LOMBWELCH\\"

path = "D:ORGANIZED\\PERIODOGRAMS_TEST_FROMFILE\\"
notes = "C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\NOTES_TEST.txt"

#path = "C:\\Users\\alzfranc\\Desktop\\ORGANIZED\\PERIODOGRAMS_FFTIG_ORIG_BACK\\"
#path = "C:\\Users\\alzfranc\\Desktop\\ORGANIZED\\PERIODOGRAMS_FFTIG_BACK\\"
#notes = "C:\\Users\\alzfranc\\Desktop\\v2\\NOTES_BACKWARD.txt"

files = os.listdir(path)[0:150]
#b = [x[:x.index('-')] for x in a]
freqs = np.linspace(55./100000,55,100000) # for periodogram of length 100k
fft_freqs = np.linspace(0, 110, 20001) # for fft

#grams  = fillPeriodogramsDictStatus(path, files, notes)
#print grams['18168']['s6']['tp9']['status']

#plotChannel('af7', grams, fft_freqs) #af7: 56

#plotAverageChannel('af7', grams, fft_freqs)
#plotAverageChannel('af7', grams, freqs)
#plotOneByOne(grams,freqs)
#plotFile("C:\\Users\\alzfranc\\Desktop\\ORGANIZED\\PERIODOGRAMS_FFTIG_ORIG_BACK\\20703-s6.txt")
#plotFile("C:\\Users\\Adam Francey\\Desktop\\LOMB EEG\\PERIODOGRAMS_BACKWARD\\201-s4.txt")
#plotFile("D:ORGANIZED\\PERIODOGRAMS_TEST_FROMFILE\\1111111111-s1_lomb_1000.txt")
#plotFile("D:ORGANIZED\\PERIODOGRAMS_TEST_FROMFILE\\1111111111-s1_lomb_orig_20amplitude.txt")
#plotFile("D:ORGANIZED\\PERIODOGRAMS_TEST_FROMFILE\\1111111111-s1.txt")
#plotSignalFile("D:\\ORGANIZED\\data\\TEST_FROMFILE\\exported - post rejections\\1111111111-s1_rejected.txt")
#plotSignalFile("D:\\ORGANIZED\\data\\TEST_FROMFILE\\exported - filtered\\1111111111-s1_filtered.txt")

# SAVE FILES
pathways = ['222222222-s1.txt','222222222-s2.txt','222222222-s3.txt','222222222-s4.txt', '222222222-s5.txt']
grams = fillPeriodogramsDictStatus(path,pathways,notes)
#grams = fillPeriodogramsDictStatus(path,files,notes)
freqs = np.linspace(float(220)/(4*100000),220,4*100000) # for periodogram of length 100k
freqs = freqs[0:int(len(freqs)/4.)]
save_band_powers(freqs, [[0,15],[15,25]], ['b1','b2'], grams, amp_spec_dens = False)



plot_tests = False

if plot_tests:
    nout = 4*100000 # was way too small at 1000!!
    Fs = 220 # sampling rate
    freqs = np.linspace(float(Fs)/nout, Fs, nout)
    freqs = freqs[0:int(len(freqs)/4.)] # only looking at frequencies 0Hz to 55Hz
    N_freqs = len(freqs)
    maxFreq = freqs[N_freqs-1]
    # fake frequency bands
    b1 = (15/maxFreq)*N_freqs
    b2 = (25/maxFreq)*N_freqs
    freqBins = [0,int(b1),int(b2)]

    def get_bands_custom(pgrams, freqBins):
        # out[i][j] = band power at channel i band j
        out = [[0,0], #tp9
               [0,0], #tp10
               [0,0], #af7
               [0,0]] #af8

        relative = [[0,0],
               [0,0],
               [0,0],
               [0,0]]

        for i in range(4):
            pgram = pgrams[i]
            for j in range(2):
                
                band = pgram[freqBins[j]:freqBins[j+1]]
                bandsum = sum(band)
                out[i][j] = bandsum

            for j in range(2):
                relative[i][j] = float(out[i][j])/sum(out[i])

        return out, relative
        


    pathways = ['222222222-s1_sine_lomb.txt','222222222-s1_sine_lombwelch.txt','222222222-s1_sine_fftignorant.txt',
                '222222222-s1_sine_fftignorant_original.txt','222222222-s1_sine_lomb_on_original.txt', '222222222-s1_sine_lomb_on_original_hamming.txt',
                '222222222-s1_sine_lombwelch_nohamming.txt']

    fig_names = ['Lomb-Scargle Periodogram', 'Lomb-Scargle-Welch Periodogram', 'Welch Periodogram',
                 'Welch Periodogram (Full Signal)', 'Lomb-Scargle Periodogram (Full Signal)', 'Lomb-Scargle Periodogram Hamming (Full Signal)',
                 'Lombwelchnohamming']


    pathways = ['222222222-s1.txt']
    fig_names = ['test']

    testpath = "D:ORGANIZED\\PERIODOGRAMS_TEST_FROMFILE\\"

    for p,n in zip(pathways, fig_names):

        pgrams = getPeriodogramFromFile(testpath+p)
        bands, relative = get_bands_custom(pgrams, freqBins)

        plt.clf()
        
        fig = plt.figure()
        fig.suptitle(n + ": b1 = %.4f, b2 = %.4f" % (relative[0][0], relative[0][1]), fontsize=12)
        plt.xlabel("Frequency (Hertz)")
        plt.ylabel("Amplitude (Volts)")
        #plotFile(testpath + p)
        plotChannelFromFile(testpath + p, 'af7')
        plt.savefig('sample pics\\' + n)
    
    
