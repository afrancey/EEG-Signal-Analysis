# use butterworth filter to filter signal noise
import scipy.signal as scisig
import os

# filter order: https://dsp.stackexchange.com/questions/37646/filter-order-rule-of-thumb

inputdir = "C:/Users/alzfr/Desktop/testeeg"
outputdir = "C:/Users/alzfr/Desktop/testeegout"
file_list = os.listdir(inputdir)

sos = scisig.butter(10, [1,30], 'bp', fs=220, output='sos')
header = "\n"

for filename in file_list:

    if "EEG" in filename:
        chandata = [[],[],[],[]]

        # open file and extract data
        with open(inputdir + "/" + filename, "r") as f:
            lines = f.readlines()
            header = lines[0]
            for line in lines[1:]:
                splitline = line.split("\t")
                for ch in range(0,4):
                    chandata[ch].append(float(splitline[ch+1]))
                    
        # filter all channels
        filtered = [scisig.sosfilt(sos, chan) for chan in chandata]

        # write to file
        with open(outputdir + "/" + filename,"w") as f:
            f.write(header)
            for sample in range(0, len(filtered[0])):
                f.write(str(sample) + "\t" +
                        str(filtered[0][sample]) + "\t" +
                        str(filtered[1][sample]) + "\t" +
                        str(filtered[2][sample]) + "\t" +
                        str(filtered[3][sample]) + "\n")
            