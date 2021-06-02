import random
import math
from numpy.lib.function_base import average

from sk_dsp_comm.fec_conv import fec_conv

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def box_muller() -> float:
       u1 = np.random.uniform()
       u2 = np.random.uniform()

       return np.sqrt(-2*np.log(u1))*np.sin(2*np.pi*u2)



if __name__ == '__main__':

       N_BITS = 10000


       snr_in_db_range = np.arange(0, 12, 0.5)

       bers = []

       
       for index in range(len(snr_in_db_range)):
              local_bers = []

              snr_in_db = snr_in_db_range[index]

              snr_linear = 10**(snr_in_db/10)

              noise_spectral_dens = 1/snr_linear

              for i in range(100):
                     

                     msg = np.random.randint(0, 2, N_BITS, int)
                     signal = 1 - 2*msg          

                     print('n0', noise_spectral_dens)

                     randoms = [box_muller() for _ in range(N_BITS)]

                     noise = math.sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                     received = signal + noise

                     demodulated = [1 if r < 0 else 0 for r in received]

                     n_errors = sum(demodulated != msg)
                     ber = n_errors/N_BITS
                     local_bers.append(ber)


              bers.append(average(local_bers))

       plt.title('BER in uncoded BPSK modulated AWGN channel')
       plt.xlabel(r'$E_b/N_0$ (dB)')
       plt.ylabel('BER')
       plt.yscale('log')
       plt.plot(snr_in_db_range, np.asarray(bers))
       plt.show()






       


