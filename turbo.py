from functools import reduce
from convolutional import ConvCode, chunks
import os
from typing import Callable, Dict, List, Optional
from enum import Enum
from itertools import product
import numpy as np
from numpy.core.numeric import Infinity
from numpy.lib.function_base import append
from numpy.lib.function_base import average
import time
from box_muller import box_muller
from math import sqrt
import csv
from interleaver import Interleaver
from rsc import RSC, modulate, pairwise_compare
import sys
import json

def flatten(list_like):
    if isinstance(list_like, (list, tuple, set, range)):
        for sub in list_like:
            yield from flatten(sub)
    else:
        yield list_like


class TurboCode():
    code: RSC
    interleaver: Interleaver

    def __init__(self, code, interleaver):
        self.code = code
        self.interleaver = interleaver

    def encode(self, sequence: List[bool]) -> List[bool]:
        output = self.code.encode(sequence)
        sys = output[::2]
        parity1 = output[1::2]

        interleaved = self.interleaver.interleave(sequence)
        parity2 = self.code.encode(interleaved)[1::2]
        # print('O2: ', output2)
        final = []

        for index, s in enumerate(sys):
            final.append(s)
            final.append(parity1[index])
            final.append(parity2[index]) 
        # print(final)
        return final

    def encode_puncture(self, sequence: List[bool]) -> List[bool]:
        encoded = self.encode(sequence)
        systematic = encoded[::3]
        parity1 = encoded[1::3]
        parity2 = encoded[2::3]

        output = []

        for index, s in enumerate(systematic):
            output.append(s)
            if index % 2 == 0:
                output.append(parity1[index])
            else:
                output.append(parity2[index])

        return output


    def decode(self, sequence: List[float], iterations):
        systematic = sequence[0::3]
        parity1 = sequence[1::3]
        parity2 = sequence[2::3]
        systematic_interleaved = self.interleaver.interleave(systematic)

        last_decoded = []
        arr_1 = []
        arr_2 = []
        arr_1.append(self.code.bcjr_decode(systematic, parity1, Le=None)[1])
        result2 = self.code.bcjr_decode(systematic_interleaved, parity2, Le=self.interleaver.interleave(arr_1[0]))
        arr_2.append(self.interleaver.deinterleave(result2[1]))
        last_decoded = self.interleaver.deinterleave(result2[0])

        for t in range(0, iterations-1):

            arr_1.append(self.code.bcjr_decode(systematic, parity1, Le=arr_2[t])[1])
            result2 = self.code.bcjr_decode(systematic_interleaved, parity2, Le=self.interleaver.interleave(arr_1[t]))
            #print('Le of 2: ', result2[1])
            arr_2.append(self.interleaver.deinterleave(result2[1]))
            last_decoded = self.interleaver.deinterleave(result2[0])
        # print('First: ', arr_1)
        # print('Second: ', arr_2)
        # for i in range(len(arr_1[0])-1):
        #     print(arr_1[0][i]-arr_1[1][i])

        return last_decoded

    def decode_punctured(self, sequence: List[float], iterations):
        systematic = sequence[0::2]
        parity1 = sequence[1::2][::2]
        parity2 = sequence[1::2][1::2]
        # print(parity1)
        # print(parity2)
        systematic_interleaved = self.interleaver.interleave(systematic)
        
        last_decoded = []
        # Filling in punctured spots with 0 values
        parity1 = list(flatten(list(zip(parity1, [0]*len(parity1)))))
        parity2 = list(flatten(list(zip([0]*len(parity2), parity2))))

        # print('p1: ', parity1)
        # print('p2: ', parity2)
        arr_1 = []
        arr_2 = []
        arr_1.append(self.code.bcjr_decode(systematic, parity1, Le=None)[1])
        result2 = self.code.bcjr_decode(systematic_interleaved, parity2, Le=self.interleaver.interleave(arr_1[0]))
        arr_2.append(self.interleaver.deinterleave(result2[1]))
        last_decoded = self.interleaver.deinterleave(result2[0])

        for t in range(0, iterations-1):

            arr_1.append(self.code.bcjr_decode(systematic, parity1, Le=arr_2[t])[1])
            result2 = self.code.bcjr_decode(systematic_interleaved, parity2, Le=self.interleaver.interleave(arr_1[t]))
            # print('Le of 2: ', result2[1])
            arr_2.append(self.interleaver.deinterleave(result2[1]))
            last_decoded = self.interleaver.deinterleave(result2[0])
        # print('First: ', arr_1)
        # print('Second: ', arr_2)
        # for i in range(len(arr_1[0])-1):
        #     print(arr_1[0][i]-arr_1[1][i])
        print(last_decoded)
        return last_decoded

def gates(register: List[bool]):
    return [register[-3] ^ register[-1]]

if __name__ == '__main__' and sys.argv[-1] == 'test':
    # random_bits = np.random.randint(0, 2, 10, bool)
    # print(random_bits)

    # def gates(input: List[bool]):
    #     return [input[0] ^ input[2]]

    # rsc = RSC(gates, [1], 3, 2)

    # coded = rsc.encode(random_bits)
    # decoded = rsc.decode(coded)
    # print(coded)
    # print(decoded)

    rsc = RSC(gates, [-2, -1], 2, 2)
    #rsc2 = RSC(gates, [1, 2], 3, 2)
    interleaver = Interleaver(int(sys.argv[1]))

    code = TurboCode(rsc, interleaver)
    random_bits = np.random.randint(0,2,int(sys.argv[1]),bool)

    

    coded = code.encode_puncture(random_bits)

    randoms = [box_muller() for _ in range(len(coded))]

    snr_linear = 10**(0/10)
    noise_variance = 0.5/snr_linear

    code.code.set_snr(snr_linear)


    noise = sqrt(noise_variance)*np.asarray(randoms)

    print(len(coded))
    modulated = modulate(coded)
    received = modulated + noise

    decoded = code.decode_punctured(received, int(sys.argv[2]))
    print(random_bits)
    print(decoded)
    n_errors = sum(1 if a != b else 0 for a, b in zip(decoded, random_bits))
    # print(coded)
    # print(decoded)
    print('Correctness: ', pairwise_compare(random_bits, decoded))
    print('N Errors: ', n_errors)

if __name__ == '__main__' and sys.argv[-1] == 'simulation':
    simulation_data = {}
    start = 0
    end = 0

    bers = []

    

    BITS = int(sys.argv[1]) # 10000

    rsc = RSC(gates, [-2, -1], 2, 2)
    interleaver = Interleaver(BITS)
    turbocode = TurboCode(rsc, interleaver)

    snr_in_db_range = np.arange(0, 8 if sys.argv[-3] != 'reduced' else 1, 0.5)

    start = time.time()

    for index in range(len(snr_in_db_range)):
        local_bers = []

        snr_in_db = snr_in_db_range[index]

        simulation_data[snr_in_db] = {}

        snr_linear = 10**(snr_in_db/10)


        noise_spectral_dens = 1/snr_linear
        noise_variance = 0.5/snr_linear

        turbocode.code.set_snr(snr_linear)

        print('n0', noise_variance)

        for i in range(int(sys.argv[2])): # 10):
            simulation_data[snr_in_db][i] = {}
            
            uncoded = np.random.randint(0, 2, BITS, bool)

            if sys.argv[-2] == 'punctured':
                msg = turbocode.encode_puncture(uncoded)
            else:
                msg = turbocode.encode(uncoded)

            signal = modulate(msg)


            randoms = [box_muller() for _ in range(len(msg))]

            noise = sqrt(noise_variance)*np.asarray(randoms)

            received = signal + noise

            if sys.argv[-2] == 'punctured':
                decoded = turbocode.decode_punctured(received, int(sys.argv[3]))
            else:
                decoded = turbocode.decode(received, int(sys.argv[3]))

            n_errors = sum(1 if a != b else 0 for a, b in zip(decoded, uncoded))
            ber = n_errors/BITS

            simulation_data[snr_in_db][i]['number_of_errors'] = int(n_errors)
            simulation_data[snr_in_db][i]['bit_error_rate'] = float(ber)

            local_bers.append(ber)

        bers.append(average(local_bers))

    end = time.time()
    print('Elapsed: ', end - start)

    print(bers)

    if not os.path.exists('turbo'):
        os.makedirs('turbo')
    if sys.argv[-2] == 'punctured':
        with open(f'turbo/punctured_turbo_iter_{sys.argv[3]}_aggregated_bers_bits_{sys.argv[1]}_n_tests_{sys.argv[2]}_{time.strftime("%m-%d-%Y-%H-%M-%S", time.localtime(end))}.csv', mode='w') as file:
            writer = csv.writer(file)
            writer.writerows(map(lambda x: [x], bers))

        with open(f'turbo/punctured_turbo_iter_{sys.argv[3]}_raw_data_bits_{sys.argv[1]}_n_tests_{sys.argv[2]}_{time.strftime("%m-%d-%Y-%H-%M-%S", time.localtime(end))}.json', mode='w') as file:
            json.dump(simulation_data, file)
    else:
        with open(f'turbo/turbo_iter_{sys.argv[3]}_aggregated_bers_bits_{sys.argv[1]}_n_tests_{sys.argv[2]}_{time.strftime("%m-%d-%Y-%H-%M-%S", time.localtime(end))}.csv', mode='w') as file:
            writer = csv.writer(file)
            writer.writerows(map(lambda x: [x], bers))

        with open(f'turbo/turbo_iter_{sys.argv[3]}_raw_data_bits_{sys.argv[1]}_n_tests_{sys.argv[2]}_{time.strftime("%m-%d-%Y-%H-%M-%S", time.localtime(end))}.json', mode='w') as file:
            json.dump(simulation_data, file)    