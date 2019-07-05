# processes EEG files

import scipy.signal as signal
import numpy as np
import os.path

class EEGSet():

        def __init__(self, originalFilename, rejectFilename, eventsFilename):

                self.okayToProcess = False
                self.error = 'None'

                self.Fs = 220 #sampling rate of our EEG data
                
                self.originalSet = self.importEEGSet(originalFilename)
                self.rejectSet = self.importEEGSet(rejectFilename)
                if self.rejectSet != "file does not exist" and self.originalSet != "file does not exist":
                        print "FILE EXISTS"
                        self.events = self.importEvents(eventsFilename)
                        self.indicatorArray = self.makeIndicatorArray(self.originalSet,self.events, self.rejectSet)
                        if self.indicatorArray != False and self.TEST_indicator(self.originalSet,self.rejectSet, self.indicatorArray):
                                self.timeSteps = self.makeTimeSteps(self.Fs, self.indicatorArray)
                                self.hammingArray = self.makeHammingArray(self.indicatorArray)
                                self.okayToProcess = True
                        else:
                                self.error = 'indicator failure'
                else:
                        self.error = 'file failure'


        def process(self, freqs, method = 'lomb'):

                N_freqs = len(freqs)
                maxFreq = freqs[N_freqs-1]
                delta_i = (4/maxFreq)*N_freqs #upper bound of delta waves
                theta_i = (8/maxFreq)*N_freqs
                alpha_i = (14/maxFreq)*N_freqs
                beta_i = (30/maxFreq)*N_freqs
                gamma_i = (50/maxFreq)*N_freqs
                freqBins = [0,int(delta_i),int(theta_i),int(alpha_i),int(beta_i),int(gamma_i)]
                
                ang_freqs = 2*np.pi*freqs

                if method == 'lomb':

                    pgrams = self.getPeriodograms_SLOW(self.rejectSet, self.timeSteps, self.hammingArray, True, ang_freqs)
                    bandpowers, relative = self.get_bands(pgrams, freqBins)

                elif method == 'lomb_on_original':

                    fullHamming = signal.hamming(len(self.originalSet[0]), sym=True)
                    #fullHamming = np.ones(len(self.originalSet[0]))
                    indicatorArray = np.ones(len(self.originalSet[0]))
                    timesteps = self.makeTimeSteps(self.Fs, indicatorArray)

                    pgrams = self.getPeriodograms_SLOW(self.originalSet, timesteps, fullHamming, True, ang_freqs)
                    bandpowers, relative = self.get_bands(pgrams, freqBins)
                    
                
                elif method == 'lombwelch':
                    windowLength = 220 # 1 second
                    windowOverlap = 110 # 0.5 seconds overlap

                    windowLength = len(self.originalSet[0])/2 # half of signal length
                    windowOverlap = windowLength/2 # windows overlap by half of their length
                    
                    pgrams = self.getPeriodograms_lombwelch(windowLength, windowOverlap, self.originalSet, self.indicatorArray, ang_freqs)
                    bandpowers, relative = 0,0
                elif method == 'fftignorant':
                    (fft_freqs, pgrams) = self.getPeriodograms_fftignorant(self.rejectSet, len(freqs))
                    bandpowers, relative = 0,0
                elif method == 'fftignorant_original':
                    (fft_freqs, pgrams) = self.getPeriodograms_fftignorant(self.originalSet,len(freqs))
                    bandpowers, relative = 0,0

                return pgrams, bandpowers, relative

        def importEEGSet(self, EEGSetFilename):

            print(os.path.isfile(EEGSetFilename))

            str_chandata = [[],[],[],[]]

            if os.path.isfile(EEGSetFilename):

                # first line is "tp9, tp10, af7, af8"
                # next n lines are the values, in the same order
                try:

                    with open(EEGSetFilename,'r') as f:
                        lines = f.readlines()

                        for line in lines[1:]:
                            line = line.strip().split('\t') #index, tp9, tp10, af7, af8

                            for i in range(4):
                                str_chandata[i].append(float(line[i+1]))

                            tp9.append(float(line[0]))
                            tp10.append(float(line[1]))
                            af7.append(float(line[2]))
                            af8.append(float(line[3]))
                    return(str_chandata)
                except ValueError:
                    print "ValueError"
                    return "file does not exist" # actually does exist but not right format (probably empty)
            else:
                    return "file does not exist"


        def importBoundaries(self, sample_boundaries_filename):

            # eventset is set of boundaries
            # latency - position of boundary
            # duration - how many samples were rejected

            # file may just contain 'none'

            sample_boundaries = []

            if os.path.isfile(sample_boundaries_filename):
                with open(sample_boundaries_filename, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        splitline = line.split(",")
                        if splitline[0] in sample_boundaries_filename:
                            sample_boundaries = [int(x) for x in splitline[1:]]


                return(sample_boundaries)
            else:
                return "file does not exist"



        def makeTimeSteps(self, Fs, indicatorArray):

            h = 1./Fs

            time = 0
            times = []

            for indicator in indicatorArray:
                time += h
                if indicator == 1:
                        times.append(time)

            return times
        
        def makeIndicatorArray(self, sample_bounds, size):
            indicatorArray = []
            for d in range(size):
                rejected = False
                for bound in range(0,len(sample_bounds),2):
                    start = min(sample_bounds[bound], sample_bounds[bound+1])
                    end = max(sample_bounds[bound], sample_bounds[bound+1])
                    if d >= start and d <=end:
                        rejected = True

                if rejected:
                    indicatorArray.append(0)
                else:
                    indicatorArray.append(1)
                    
            return indicatorArray

        def makeHammingArray(self, indicatorArray):

                numSamples = len(indicatorArray)
                fullHamming = signal.hamming(numSamples, sym=True)

                trimmedHamming = []
                for sampleNum in range(len(indicatorArray)):
                        indicator = indicatorArray[sampleNum]
                        if indicator == 1:
                                trimmedHamming.append(fullHamming[sampleNum])

                return trimmedHamming

        def getPeriodograms_lombwelch(self, windowLength, windowOverlap, original, indicatorArray, ang_freqs, normalize = True):

            N = len(indicatorArray)

            pgrams = []
            for ch in original:

                ch_array = np.array(ch)

                pgramSum = np.zeros(len(ang_freqs))
                pgramCount = 0

                for i in range(0,N, windowOverlap):
                    if windowLength + i > N:
                        # this should be last window
                        endIndex = N
                    else:
                        endIndex = i+windowLength
                        
                    indicator = indicatorArray[i:endIndex]
                    if sum(indicator) > 0:
                        # if there are any good samples
                        samples = ch_array[i:endIndex]
                        
                        t = np.array(self.makeTimeSteps(220, indicator))
                        ham = np.array(self.makeHammingArray(indicator))

                        trimmedSamples = np.array([])
                        count = 0
                        for ind in indicator:
                            if ind:
                                trimmedSamples = np.append(trimmedSamples,samples[count])
                            count += 1

                        y = trimmedSamples*ham
                        y = y - np.mean(y)

                        pgram = signal.lombscargle(t,y,ang_freqs)
                        if normalize:
                            pgram = np.sqrt(4*(pgram/len(pgram)))
                        pgramSum = pgramSum + pgram
                        pgramCount += 1
                   #else: no samples within this time window, do nothing
                avgPgram = pgramSum/pgramCount
                
                #deprecated: now normalizing intermediate periodograms
                #if normalize:
                #    avgPgram = np.sqrt(4*(avgPgram/len(avgPgram)))
                    
                pgrams.append(avgPgram)
            return pgrams

        def get_bands(self, pgrams, freqBins):

            # out[i][j] = band power at channel i band j
            out = [[0,0,0,0,0], #tp9
                   [0,0,0,0,0], #tp10
                   [0,0,0,0,0], #af7
                   [0,0,0,0,0]] #af8

            relative = [[0,0,0,0,0],
                   [0,0,0,0,0],
                   [0,0,0,0,0],
                   [0,0,0,0,0]]

            for i in range(4):
                pgram = pgrams[i]
                for j in range(5):
                    
                    band = pgram[freqBins[j]:freqBins[j+1]]
                    bandsum = sum(band)
                    out[i][j] = bandsum

                for j in range(5):
                    relative[i][j] = float(out[i][j])/sum(out[i])

            return out, relative
