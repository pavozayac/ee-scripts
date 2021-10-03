from math import sqrt
from random import uniform
from numpy.lib.function_base import average
from box_muller import box_muller

# import matplotlib.pyplot as plt
import numpy as np

import sys
import csv
from utils import chunks
import os 
import json 
import time

# plt.xlabel(r'$E_b/N_0$ (dB)')
# plt.ylabel('BER')
# plt.yscale('log')

# all_ber_arrays = []


#
# Uncoded BPSK
#

if __name__ == '__main__':  
              
       N_BITS = 10000 * 100

       simulation_data = {}

       snr_in_db_range = np.arange(0, 9, 0.5)

       uncoded_bers = []

       start = time.time()

       for index in range(len(snr_in_db_range)):
              local_bers = []

              snr_in_db = snr_in_db_range[index]

              simulation_data[snr_in_db] = {}

              snr_linear = 10**(snr_in_db/10)

              noise_spectral_dens = 1/snr_linear

              for i in range(10):
                     simulation_data[snr_in_db][i] = {}

                     msg = np.random.randint(0, 2, N_BITS, int)

                     signal = 1 - 2*msg          

                     print('n0', noise_spectral_dens)

                     randoms = [box_muller() for _ in range(N_BITS)]

                     noise = sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                     received = signal + noise

                     demodulated = [1 if r < 0 else 0 for r in received]

                     n_errors = sum(demodulated != msg)
                     ber = n_errors/N_BITS

                     simulation_data[snr_in_db][i]['number_of_errors'] = int(n_errors)
                     simulation_data[snr_in_db][i]['bit_error_rate'] = float(ber)   

                     local_bers.append(ber)

              uncoded_bers.append(average(local_bers))

       end = time.time()
       print('Elapsed: ', end - start)

       if not os.path.exists('channel'):
              os.makedirs('channel')

       with open(f'channel/channel_aggregated_bers_{end}.csv', mode='w') as file:
              writer = csv.writer(file)
              writer.writerows(map(lambda x: [x], uncoded_bers))

       with open(f'channel/channel_raw_data_{end}.json', mode='w') as file:
              json.dump(simulation_data, file)