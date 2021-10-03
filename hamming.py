from typing import List
from functools import lru_cache

@lru_cache(maxsize=None)
def encode_block(n: tuple):
    m = list(n)
    m.insert(0, 0)
    m.insert(1, 0)
    m.insert(3, 0)
    m[0] = (m[0]+m[2]+m[4]+m[6])%2
    m[1] = (m[1]+m[2]+m[5]+m[6])%2
    m[3] = (sum(m[3:7]))%2

    return m
 
def hamming_encode (arr: list):
    if len(arr) % 4 != 0:
        raise Exception("Message not aligned to 4 bits")
    
    blocks = []
    for i in range(int(len(arr)/4)):
        blocks.append(tuple(arr[i*4:i*4+4]))

    return_blocks = []

    for m in blocks:
        return_blocks.append(encode_block(m))
        
    return [bit for block in return_blocks for bit in block]

@lru_cache(maxsize=None)
def correct_block(n: tuple):
    m = list(n)
    syndrome = []
    syndrome.append((sum(m[3:7]))%2)
    syndrome.append((m[1]+m[2]+m[5]+m[6])%2)
    syndrome.append((m[0]+m[2]+m[4]+m[6])%2)
    if syndrome != [0,0,0]:
        position = -1
        for bit in syndrome:
            position = (position<<1) | bit
        if m[position] == 0:
            m[position] = 1
        else:
            m[position] = 0
            
    return m

def hamming_correct(arr: list):
    blocks = []
    for i in range(int(len(arr)/7)):
        blocks.append(tuple(arr[i*7:i*7+7]))

    return_blocks = []

    for m in blocks:
        return_blocks.append(correct_block(m))
            
    return [bit for block in return_blocks for bit in block]

def hamming_decode(arr: list):
    
    blocks = []
    for i in range(int(len(arr)/7)):
        blocks.append(arr[i*7:i*7+7])

    output = []
    for m in blocks:
        output.append([m[2]]+m[4:7])

    return [bit for block in output for bit in block]


from utils import chunks
from box_muller import box_muller
import sys
import time
import csv
import json
import os
import numpy as np

if __name__ == '__main__':

    HAMMING_BITS = 10000 * 100
    snr_in_db_range = np.arange(0, 9, 0.5)
    hamming_bers = []

    simulation_data = {}


    start = time.time()

    for index in range(len(snr_in_db_range)):
        local_bers = []

        snr_in_db = snr_in_db_range[index]

        simulation_data[snr_in_db] = {}

        snr_linear = 10**(snr_in_db/10)

        noise_spectral_dens = 1/snr_linear

        for i in range(10):
            simulation_data[snr_in_db][i] = {}
            
            uncoded = np.random.randint(0, 2, HAMMING_BITS, int)

            #msg = hamming.hamm_encoder(uncoded)
            msg = hamming_encode(uncoded)

            signal = 1 - 2*np.array(msg)

            print('n0', noise_spectral_dens)

            randoms = [box_muller() for _ in range(len(msg))]

            noise = np.sqrt(noise_spectral_dens/2)*np.asarray(randoms)

            received = signal + noise

            demodulated = [1 if r < 0 else 0 for r in received]

            #decoded = hamming.hamm_decoder(np.asarray(demodulated).astype(int))
            corrected = hamming_correct(demodulated)
            decoded = hamming_decode(corrected)

            n_errors = sum(1 if a != b else 0 for a, b in zip(decoded, uncoded))
            ber = n_errors/HAMMING_BITS

            simulation_data[snr_in_db][i]['number_of_errors'] = int(n_errors)
            simulation_data[snr_in_db][i]['bit_error_rate'] = float(ber)

            local_bers.append(ber)

        hamming_bers.append(np.average(local_bers))

    print(hamming_bers)
    end = time.time()
    print('Elapsed: ', end - start)


    if not os.path.exists('hamming'):
        os.makedirs('hamming')

    with open(f'hamming/hamming_aggregated_bers_{end}.csv', mode='w') as file:
            writer = csv.writer(file)
            writer.writerows(map(lambda x: [x], hamming_bers))

    with open(f'hamming/hamming_raw_data_{end}.json', mode='w') as file:
        json.dump(simulation_data, file)