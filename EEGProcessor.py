# processes EEG files

import scipy.signal as signal
import numpy as np
import os.path

class EEGSet():

    def __init__(self, originalFilename, eventsFilename):

        # Chooseable Parameters:
        
        # UPPER BOUNDS of each band
        HZ_DELTA = 4
        HZ_THETA = 8 
        HZ_ALPHA = 14
        HZ_BETA = 30
        self.HZ_ALL_BANDS = [HZ_DELTA, HZ_THETA, HZ_ALPHA, HZ_BETA]

        self.Fs = 220 #sampling rate of our EEG data
        self.windowLength = 220 # 1 second
        self.windowOverlap = 110 # 0.5 seconds overlap

        self.samplingTime = 60 # seconds
        self.num_samples_full = self.samplingTime*self.Fs

        # set the frequencies to evaluate
        # range: because of Nyquist limit, we may sample frequencies up to Fs/2 Hz
        # grid size: 1/(N * 1/Fs) = Fs/N
        # freqs = [0, Fs/N, 2Fs/N, 3Fs/N,..., (N-1)Fs/N]
        # also, ang_freqs = 2*np.pi*freqs
        # See: 
        # If
        # N = num samples,
        # df = spacing between frequencies in periodogram
        # dt = time between samples
        # then we require
        # df >= 1/(N*dt)
        # source: Jacob T. VanderPlas, "Understanding the Lomb-Scargle Periodogram", pg 14
        #           https://arxiv.org/pdf/1703.09824.pdf
        self.freqs = np.linspace(0, (self.num_samples_full - 1)*self.Fs/self.num_samples_full, self.num_samples_full)

        # cut off freqs at 30 Hz
        # index at which we hit 30 Hz or more:
        # freqs[i] = i*Fs/N, 0 <= i <= N-1
        # we wish for freqs[i] >= 30 Hz
        # -> i*Fs/N >= 30
        # -> i >= 30*N/Fs
        self.freqs = self.freqs[:int(30*self.num_samples_full/self.Fs)]


        self.windowTimesteps = np.linspace(0, 1-1/self.windowLength, self.windowLength)
        self.windowHamming = signal.hamming(self.windowLength, sym=True)


        # Begin init
                        
        self.okayToProcess = False
        self.error = 'None'

        self.filename = originalFilename
        
        self.originalSet = self.importEEGSet(originalFilename)
        if self.originalSet != "file does not exist":
                print("FILE EXISTS")
                self.sample_boundaries = self.importBoundaries(eventsFilename)
                self.indicatorArrays = [self.makeIndicatorArray(self.sample_boundaries[ch], len(self.originalSet[ch])) for ch in range(4)]
                self.okayToProcess = True
        else:
                self.error = 'file failure'


    def process(self):

            # calculate periodograms of each channel
            pgrams = self.getPeriodograms_lombwelch()

            # sum values in each frequency bin
            # first change HZ_ALL_BANDS into list of freq bin boundary array indices
            N_freqs = len(self.freqs)
            maxFreq = self.freqs[N_freqs-1]
            freqBins = [0] + [int((FREQ/maxFreq)*N_freqs) for FREQ in self.HZ_ALL_BANDS]

            # get powers between each index
            bandpowers, relative = self.get_bands(pgrams, freqBins)

            return pgrams, bandpowers, relative

    def importEEGSet(self, EEGSetFilename):

        print(os.path.isfile(EEGSetFilename))

        str_chandata = [[],[],[],[]]

        if os.path.isfile(EEGSetFilename):

            # first line is "index,tp9, tp10, af7, af8"
            # next n lines are the values, in the same order
            try:

                with open(EEGSetFilename,'r') as f:
                    lines = f.readlines()

                    for line in lines[1:]:
                        line = line.strip().split('\t') #index, tp9, tp10, af7, af8

                        for i in range(4):
                            str_chandata[i].append(float(line[i+1]))

                return(str_chandata)
            except ValueError:
                print("ValueError")
                return "file does not exist" # actually does exist but not right format (probably empty)
        else:
                return "file does not exist"


    def importBoundaries(self, sample_boundaries_filename):
        # sample boundaries file has lines that look like this:
        
        # "name, # # #..., # # #..., # # #..., # # #..."
        # where commas separate channels.
        # in each channel, boundaries come in pairs

        sample_boundaries = []

        if os.path.isfile(sample_boundaries_filename):
            with open(sample_boundaries_filename, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    splitline = line.split(",")
                    if splitline[0] in self.filename:
                        # we have found our boundaries in the boundary file
                        for chan in splitline[2:]: #first two items in .csv row are taken up by filename in this case (remember extra comma)
                            boundaries_int = [int(x) for x in chan.split(" ")[:-1]] # last character is space, not an int for boundary marker
                            sample_boundaries.append(boundaries_int)


            return(sample_boundaries)
        else:
            return "file does not exist"


    def trim(self, series, indicatorArray):
        # trims a list (series) according to indicatorArray
        # series: time series whose interval matches up with indicatorArray

        trimmedSeries = []
        for count in len(indicatorArray):
            if indicatorArray[count] == 1:
                trimmedSeries.append(series[count])

        return(trimmedSeries)

    def trim_nparray(self, nparray, indicatorArray):
        # trims a list (series) according to indicatorArray
        # series: time series whose interval matches up with indicatorArray

        trimmedSeries = np.array([])
        for count in range(len(indicatorArray)):
            if indicatorArray[count] == 1:
                np.append(trimmedSeries,nparray[count])

        return(trimmedSeries)
    
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

    def lombscarglewelch(self, input_series, input_indicatorArray):

        # Calculates the average periodogram over a sequence on overlapping time-windows

        # truncate the series and indicatorArray to nearest multiple of windowOverlap
        N = int(len(input_series)/self.windowOverlap)*self.windowOverlap       
        series = input_series[0:N - 1]
        indicatorArray = input_indicatorArray[0:N - 1]

        # we will be element-wise summing periodograms as we go and taking count
        pgramSum = np.zeros(len(self.freqs))
        pgramCount = 0

        # get periodogram for each window, counting the ending of each window
        for i in range(self.windowLength,N, self.windowOverlap):

            # determine sample interval
            startIndex = i - self.windowLength
            endIndex = i - 1

            # get indicator array for this interval
            indicator = indicatorArray[startIndex:endIndex]

            # check if there are any good samples
            if sum(indicator) > 0:

                # get interval of samples                      
                samples = np.array(series[startIndex:endIndex])

                # trim the relevant arrays
                t = self.trim_nparray(self.windowTimesteps,indicator)
                ham = self.trim_nparray(self.windowHamming,indicator)
                trimmedSamples = self.trim_nparray(samples,indicator)

                # apply Hamming window element-wise
                y = trimmedSamples*ham

                # calculate periodogram

                # to do: change to using astroML
                ###########################################
                # pgram = signal.lombscargle(t,y,ang_freqs)
                # HACK
                pgram = np.zeros(len(self.freqs))
                ##########################################

                # add to running sums
                pgramSum = pgramSum + pgram
                pgramCount += 1
                                      
           #else: no samples within this time window, do nothing

        # find average periodogram
        avgPgram = pgramSum/pgramCount
        return(avgPgram)
            
    def getPeriodograms_lombwelch(self):

        pgrams = []
        for channel in range(len(self.originalSet)):
            avgPgram = self.lombscarglewelch(self.originalSet[channel],self.indicatorArrays[channel])
            pgrams.append(avgPgram)
        return pgrams

    def getPeriodograms_lombwelch_OLD(self, windowLength, windowOverlap, original, indicatorArrays, ang_freqs, normalize = True):

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
            for j in range(len(freqBins) - 1):

                band = pgram[freqBins[j]:freqBins[j+1]]
                bandsum = sum(band)
                out[i][j] = bandsum

            for j in range(5):
                relative[i][j] = float(out[i][j])/sum(out[i])

        return out, relative
