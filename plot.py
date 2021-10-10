from os import read
from numpy import arange
import matplotlib.pyplot as plot
from sys import argv
from csv import reader


filenames = argv[1:-1]
snr_in_db_range = arange(0, int(argv[-1]), 0.5)

plot.figure()
plot.subplot()
plot.xlabel(r'$E_b/N_0$ (dB)')
plot.ylabel('BER')

for filename in filenames:
    with open(filename, newline='\n') as file:
        
        plot.yscale('log')
        plot.xscale('linear')

        ber_reader = reader(file, delimiter=' ')

        bers = []

        for row in list(ber_reader)[:int(argv[-1])*2]:
            print(row)
            print(row[0])
            bers.append(float(row[0]))

        print(bers)

        plot.plot(snr_in_db_range, bers, label=filename)

plot.title('BER in the BPSK modulated AWGN channel')
plot.legend()
plot.show()