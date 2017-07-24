# gets periodograms of a file with 4 columns
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
tp9 = []
tp10 = []
af7 = []
af8 = []

#filename = 'C:\\Users\\alzfranc\\Desktop\\ORGANIZED\\data\\BACKWARD\\EEG\\split by stop\\1446314850918-20703\\1446314850918-20703-s6.txt'
#filename = 'C:\\Users\\alzfranc\\Desktop\\ORGANIZED\\data\\BACKWARD\\EEG\\split by stop\\1447527918014-201\\1447527918014-201-s4.txt'

filename = 'D:\\ORGANIZED\\data\\BACKWARD\\EEG\\split by stop\\1447527918014-201\\1447527918014-201-s4.txt'

with open(filename, 'r') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip().split('\t')
        tp9.append(float(line[0]))
        tp10.append(float(line[1]))
        af7.append(float(line[2]))
        af8.append(float(line[3]))


(freqlist, pgrams) = signal.periodogram(np.array([tp9,tp10,af7,af8]), fs = 220, window = 'hamming',
                                   nfft = 40000, detrend = 'constant', return_onesided = True,
                                        scaling = 'density', axis = -1)

fft_freqs = np.linspace(0, 110, 20001)


for i in range(1,5):
    plt.subplot(2,2,i)
    plt.plot(fft_freqs[:len(fft_freqs)/2.],pgrams[i-1][:len(fft_freqs)/2.])

plt.show()
