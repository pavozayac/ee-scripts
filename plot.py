from os import read
from numpy import arange
import matplotlib.pyplot as plot
from sys import argv
from csv import reader


filename = argv[1]
snr_in_db_range = arange(0, 9, 0.5)

with open(filename, newline='\n') as file:
    plot.figure()
    plot.subplot()
    plot.xlabel(r'$E_b/N_0$ (dB)')
    plot.ylabel('BER')
    plot.yscale('log')
    plot.xscale('linear')

    ber_reader = reader(file, delimiter=' ')

    bers = []

    for row in ber_reader:
        print(row[0])
        bers.append(float(row[0]))

    print(bers)

    plot.plot(snr_in_db_range, bers, label=filename)

    plot.title('BER in the BPSK modulated AWGN channel')
    plot.legend()
    plot.show()