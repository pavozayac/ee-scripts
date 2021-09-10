from math import sqrt
from random import uniform
from numpy.lib.function_base import average
from scipy.interpolate import BSpline, make_interp_spline
from box_muller import box_muller

import matplotlib.pyplot as plt
import numpy as np

import sys
import csv

def chunks(lst: list, n: int):
       for i in range(0, len(lst), n):
              yield lst[i:i + n]



plt.xlabel(r'$E_b/N_0$ (dB)')
plt.ylabel('BER')
plt.yscale('log')

all_ber_arrays = []


#
# Uncoded BPSK
#


N_BITS = 10000 * 100


snr_in_db_range = np.arange(0, 9, 0.5)

uncoded_bers = []

if sys.argv[1] == 'uncoded' or sys.argv[1] == 'all':

       for index in range(len(snr_in_db_range)):
              local_bers = []

              snr_in_db = snr_in_db_range[index]

              snr_linear = 10**(snr_in_db/10)

              noise_spectral_dens = 1/snr_linear

              for i in range(10):
                     

                     msg = np.random.randint(0, 2, N_BITS, int)
                     signal = 1 - 2*msg          

                     print('n0', noise_spectral_dens)

                     randoms = [box_muller() for _ in range(N_BITS)]

                     noise = sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                     received = signal + noise

                     demodulated = [1 if r < 0 else 0 for r in received]

                     n_errors = sum(demodulated != msg)
                     ber = n_errors/N_BITS
                     local_bers.append(ber)


              uncoded_bers.append(average(local_bers))

       plt.plot(snr_in_db_range, uncoded_bers, 'x-r', label = 'Uncoded BPSK')
       with open('uncoded.csv', mode='w') as file:
              writer = csv.writer(file)
              writer.writerows(map(lambda x: [x], uncoded_bers))



all_ber_arrays.append(uncoded_bers)


#
#      Hamming code
#

from sk_dsp_comm.fec_block import fec_hamming
import time

HAMMING_BITS = 10000 * 100


snr_in_db_range = np.arange(0, 9, 0.5)

hamming_bers = []

#hamming = fec_hamming(3)
from hamming import hamming_encode, hamming_correct, hamming_decode

if sys.argv[1] == 'hamming' or sys.argv[1] == 'all':
       start = time.time()

       for index in range(len(snr_in_db_range)):
              local_bers = []

              snr_in_db = snr_in_db_range[index]

              snr_linear = 10**(snr_in_db/10)

              noise_spectral_dens = 1/snr_linear

              for i in range(10):
                     
                     
                     uncoded = np.random.randint(0, 2, HAMMING_BITS, int)

                     #msg = hamming.hamm_encoder(uncoded)
                     msg = hamming_encode(uncoded)

                     signal = 1 - 2*np.array(msg)

                     print('n0', noise_spectral_dens)

                     randoms = [box_muller() for _ in range(len(msg))]

                     noise = sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                     received = signal + noise

                     demodulated = [1 if r < 0 else 0 for r in received]

                     #decoded = hamming.hamm_decoder(np.asarray(demodulated).astype(int))
                     corrected = hamming_correct(demodulated)
                     decoded = hamming_decode(corrected)

                     n_errors = sum(1 if a != b else 0 for a, b in zip(decoded, uncoded))
                     ber = n_errors/HAMMING_BITS
                     local_bers.append(ber)


              hamming_bers.append(average(local_bers))
       end = time.time()
       print('Elapsed: ', end - start)

       all_ber_arrays.append(hamming_bers)
       plt.plot(snr_in_db_range, hamming_bers, 'x-b', label='Hamming(7,4)')
       with open('hamming.csv', mode='w') as file:
              writer = csv.writer(file)
              writer.writerows(map(lambda x: [x], hamming_bers))

#
#      BCH code
#

from bch_lru import bch15_7

bch_code = bch15_7()

BCH_BITS = 9996 * 100

snr_in_db_range = np.arange(0, 9, 0.5)

bch_bers = []

if sys.argv[1] == 'bch' or sys.argv[1] == 'all':
       start = time.time()

       for index in range(len(snr_in_db_range)):
              local_bers = []

              snr_in_db = snr_in_db_range[index]

              snr_linear = 10**(snr_in_db/10)

              noise_spectral_dens = 1/snr_linear

              for i in range(10):
                     
                     
                     uncoded = np.random.randint(0, 2, BCH_BITS, int)

                     msg = []

                     for chunk in chunks(uncoded, 7):
                            msg += bch_code.encode(tuple(chunk))

                     signal = 1 - 2*np.asarray(msg)

                     print('n0', noise_spectral_dens)

                     randoms = [box_muller() for _ in range(len(msg))]

                     noise = sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                     received = signal + noise

                     demodulated = [1 if r < 0 else 0 for r in received]
                     
                     decoded = []

                     for dem_chunk in chunks(demodulated, 15):
                            decoded += bch_code.decode(tuple(bch_code.correct(tuple(dem_chunk))))

                     n_errors = sum(decoded != uncoded)
                     ber = n_errors/BCH_BITS
                     local_bers.append(ber)


              bch_bers.append(average(local_bers))
       end = time.time()
       print('Elapsed: ', end - start)

       all_ber_arrays.append(bch_bers)
       plt.plot(snr_in_db_range, bch_bers, 'x-b', label='BCH(15,7,2)')

       with open('bch.csv', mode='w') as file:
              writer = csv.writer(file)
              writer.writerows(map(lambda x: [x], bch_bers))

#
#      Reed-Solomon
#

from rs_lru import RS

rs_code = RS(7,3)

RS_BITS = 9999 * 100

snr_in_db_range = np.arange(0, 9, 0.5)

rs_bers = []

if sys.argv[1] == 'rs' or sys.argv[1] == 'all':
       start = time.time()

       for index in range(len(snr_in_db_range)):
              local_bers = []

              snr_in_db = snr_in_db_range[index]

              snr_linear = 10**(snr_in_db/10)

              noise_spectral_dens = 1/snr_linear

              for i in range(10):
                     uncoded = np.random.randint(0, 2, RS_BITS, int)

                     msg = []

                     for chunk in chunks(uncoded, 9):
                            msg += rs_code.encode(tuple(chunk))

                     signal = 1 - 2*np.asarray(msg)

                     print('n0', noise_spectral_dens)

                     randoms = [box_muller() for _ in range(len(msg))]

                     noise = sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                     received = signal + noise

                     demodulated = [1 if r < 0 else 0 for r in received]
                     
                     decoded = []

                     for dem_chunk in chunks(demodulated, 21):
                            decoded += rs_code.decode(tuple(rs_code.correct(tuple(dem_chunk))))

                     n_errors = sum(decoded != uncoded)
                     ber = n_errors/RS_BITS
                     local_bers.append(ber)

              rs_bers.append(average(local_bers))

       print(rs_bers)
       end = time.time()
       print('Elapsed: ', end - start)

       all_ber_arrays.append(rs_bers)
       plt.plot(snr_in_db_range, rs_bers, 'x-g', label='RS(7,3)')

       with open('rs.csv', mode='w') as file:
              writer = csv.writer(file)
              writer.writerows(map(lambda x: [x], rs_bers))



plt.title('BER in the BPSK modulated AWGN channel')
plt.legend()
plt.show()








