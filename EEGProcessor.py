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

            print os.path.isfile(EEGSetFilename)

            if os.path.isfile(EEGSetFilename):

                # first line is "tp9, tp10, af7, af8"
                # next n lines are the values, in the same order
                try:

                    with open(EEGSetFilename,'r') as f:
                        lines = f.readlines()

                        tp9 = []
                        tp10 = []
                        af7 = []
                        af8 = []

                        for line in lines[1:]:
                            line = line.strip().split('\t') #tp9, tp10, af7, af8

                            tp9.append(float(line[0]))
                            tp10.append(float(line[1]))
                            af7.append(float(line[2]))
                            af8.append(float(line[3]))
                    return [tp9,tp10,af7,af8]
                except ValueError:
                    print "ValueError"
                    return "file does not exist" # actually does exist but not right format (probably empty)
            else:
                    return "file does not exist"


        def importEvents(self, eventSetFilename):

            # eventset is set of boundaries
            # latency - position of boundary
            # duration - how many samples were rejected

            # file may just contain 'none'

            events = []

            if os.path.isfile(eventSetFilename):
                    print True

                    with open(eventSetFilename, 'r') as f:
                        lines = f.readlines()

                        if lines[0] != 'none': # else events stays = [] 

                            # check that all events are boundaries
                            for letter in lines[1].strip():

                                if letter != 'b':
                                    print "non-event found"
                                    return False



                            for line in lines[3:]:
                                line = line.strip().split('\t') # latentcy, duration

                                events.append(line)


                    return events
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
        
        def makeIndicatorArray(self, originalSet, eventSet, rejectSet):

                # relevent arrays and such
                # arbitrarily using tp9 since indicator array is universal over all channels
                orig = originalSet[0]
                reject = rejectSet[0]

                # first break the rejectSet into chunks
                chunks = []
                prev = 0
                for event in eventSet:
                        latency = int(float(event[0]))
                        chunks.append(reject[prev:latency])
                        prev = latency
                chunks.append(reject[prev:])

                # given the chunks, construct an array of length len(orig) such that
                # indicatorArray[i] = 1 where orig[i] is located within a chunk
                indicatorArray = [1]*len(orig)
                i = 0
                lastKnownGoodi = i
                for chunk in chunks:
                        l = len(chunk)
                        foundSubset = False
                        while not foundSubset:

                                if i >= len(indicatorArray):

                                        # this catches a chunk due to overlap bug in EEGLAB
                                        print "can't find chunk number " + str(chunks.index(chunk)) + " out of " + str(len(chunks)) #
                                        print "chunk length: " + str(len(chunk))
                                        
                                        i = lastKnownGoodi
                                        foundSubset=True
                                        return False
                                elif orig[i:i+l] == chunk:
                                        # we have found the subset corresponding to the chunk
                                        # list initialized as ones - no changes needed
                                        # skip the index forward
                                        i = i +  l
                                        lastKnownGoodi = i
                                        foundSubset = True
                                else:
                                        # we are not inside the subset
                                        # this sample is in a rejected portion
                                        # set it to 0 to indicate as such
                                        indicatorArray[i] = 0
                                        i = i + 1

                # now do something special
                # if there is a boundary at the edge of the set
                # EEGLAB will NOT record it as such
                # therefore, crop indicator accordingly based on length of rejects

                N_reject = len(reject)
                N_indicator = len(indicatorArray)
                N_ones = sum(indicatorArray)
                indicatorArray = indicatorArray[:N_indicator - (N_ones-N_reject)]



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


        def getPeriodograms_SLOW(self, rejectSet, times, window, normalize, ang_freqs):
            pgrams = []
            t = np.array(times) # times
            hammingWindow = np.array(window)
            n = len(t)

            print "n: " + str(n)

            for ch in rejectSet:


                y = np.array(ch) # samples
                y = y*hammingWindow # apply hamming window
                y = y - np.mean(y) # demeaned samples ******** IS THIS NECESSARY?? Answer: yes to avoid low frequency noise

                pgram = signal.lombscargle(t, y, ang_freqs)
                if normalize:
                    pgram = np.sqrt(4*(pgram/n))

                print "len(pgram): " + str(len(pgram))

                pgrams.append(pgram)

            return pgrams

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

        def getPeriodograms_fftignorant(self,rejectSet, N):
                #return signal.periodogram(np.array(rejectSet), fs = 220, window = 'hamming',
                #                   nfft = 40000, detrend = 'constant', return_onesided = True,
                #                        scaling = 'density', axis = -1)

                # using hamming window seems to throw things off quite a bit?
                #(frequencies, pgram) = signal.periodogram(np.array(rejectSet), fs = self.Fs, window = 'hamming',
                                  #detrend = 'constant', axis = -1)
                (frequencies, pgram) = signal.periodogram(np.array(rejectSet), fs = self.Fs, nfft = 4*N,
                                  detrend = 'constant', return_onesided = False, axis = -1)

                print pgram.shape
                return (frequencies[:N], np.sqrt(4*self.Fs*pgram/len(rejectSet[0]))[:,:N])
                #return (frequencies, np.sqrt(self.Fs*pgram/len(pgram)))
                


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

        def TEST_indicator(self, originalSet,rejectSet, indicatorArray):

                # tests the indicatorArray by trimming the original set and comparing to the rejected

                original = originalSet[0]
                rejected = rejectSet[0]

                #print indicatorArray
                


                with open('testfile2.txt','w') as f:
                        for i in range(len(indicatorArray)):
                                f.write(str(i) + '\t' + str(original[i]) + '\t' + str(indicatorArray[i]) + '\n')



                trimmedOriginal = []

                rejCount = 0

                #print zip([x for x in range(len(indicatorArray))],original)

                for sampleNum in range(len(indicatorArray)):
                        indicator = indicatorArray[sampleNum]

                        if indicator == 1:
                                
                                trimmedOriginal.append(original[sampleNum])

##                                try:
##                                        if original[sampleNum] != rejected[rejCount]:
##                                                print "MISMATCH AT SAMPLE NUM " + str(sampleNum)
##                                        rejCount += 1
##                                except IndexError:
##                                        print "INDEXEROR AT SAMPLE NUM " + str(sampleNum)

                with open("testfile.txt",'w') as f:
                        falseCount = 0
                        for i in range(len(trimmedOriginal)):
                                same = trimmedOriginal[i] == rejected[i]
                                if not same:
                                        #falseCount = falseCount + 1
                                        falseCount+=1
                                f.write(str(i) + "\t" + str(trimmedOriginal[i]) + "\t" + str(rejected[i]) + "\t" + str(same) + "\n")

                        if falseCount: # == if falseCount != 0

                                print 'falseCount: ' + str(falseCount)


                return trimmedOriginal == rejected
