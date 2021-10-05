from numpy import random
from numpy.lib.arraysetops import isin
from convolutional import ConvCode, chunks
import os
from typing import Callable, Dict, List, Optional, Tuple
from enum import Enum
from itertools import product
import numpy as np
from numpy.core.numeric import Infinity
from numpy.lib.function_base import append
from numpy.lib.function_base import average
from box_muller import box_muller
from math import sqrt
import csv
from interleaver import Interleaver
import sys
import datetime
from dataclasses import dataclass

class RSCMode(Enum):
    SYSTEMATIC = 'systematic'
    NONSYSTEMATIC = 'nonsystematic'

class RSCDecoder(Enum):
    BCJR = 'bcjr'
    VITERBI = 'viterbi'

def modulate(a: bool) -> float:
    if isinstance(a, bool):
        return 1 - 2*float(a)
    else:
        return [1-2*float(b) for b in a]

#This gamma is defined specifically for BPSK modulation
def compute_gamma(signal: float, parity: float, prototype_output: List[float], a_priori = 0.5, No = 2) -> float:
    # print(signal, parity, prototype_output)
    # return a_priori*np.exp(
    #     -1*((signal-prototype_output[0])**2 + (parity-prototype_output[1])**2)/2/variance
    # )

    return 2/No*(prototype_output[0]*signal+prototype_output[1]*parity) + np.log(a_priori)

def max_star(x, y):
    return max(x, y) + np.log(1 + np.exp(-1*np.abs(y-x)))


class RSC(ConvCode):
    # RSC should always have one gate less than its output

    def __init__(self, gates, recursive_indices: List[int], reg_size = 3, out_size = 2, mode = RSCMode.SYSTEMATIC, terminate = False, decoder = RSCDecoder.VITERBI):
        self.trellis = {}
        self.reg_size = reg_size
        self.out_size = out_size
        self.gates = gates
        self.register = [False]*reg_size
        self.recursive_indices = recursive_indices
        self.mode = mode
        self.terminate = terminate
        self.decoder = decoder

        if decoder == RSCDecoder.VITERBI:
            self.fill_trellis()
        else:
            self.fill_trellis_bcjr()

    def fill_trellis_bcjr(self):
        posregs = product([modulate(False), modulate(True)], repeat=self.reg_size)

        for reg in posregs:
            self.trellis[reg] = {}

            for next_bit in [modulate(False), modulate(True)]:
                self.register = list(reg).copy()

                output = self.push_reg(next_bit)

                self.trellis[reg][next_bit] = {
                    'next_state': tuple(self.register),
                    'output': output
                }

    def push_reg(self, bit: bool):
        recursive_bit = bit
        
        for i in self.recursive_indices:
            recursive_bit = recursive_bit != self.register[i]

        self.register = self.register[1:] + [recursive_bit]
        output = [bit] + self.gates(self.register)
        return output

    def encode(self, bits: List[bool]) -> List[bool]:
        self.register = [False]*self.reg_size

        output: List[int] = []

        for bit in bits:
            output += self.push_reg(bit)

        if self.terminate == True:
            for _ in range(self.reg_size):
                output += self.push_reg(self.register[self.recursive_indices[0]])

        # print(self.register)

        return output

    

    #takes a list of floating-point numbers and returns the log likelihood ratio, optionally takes a list of previous LLRs
    def bcjr_decode(self, systematic_sequence: List[float], parity_sequence: List[float], llrs: Optional[List[float]] = None, variance = 1) -> List[float]:

        @dataclass
        class TimeslotState():
            alpha: float
            beta: float
            gamma: Dict[bool, float]

            def __repr__(self) -> str:
                return str({
                    'alpha': self.alpha,
                    'beta': self.beta,
                    'gamma': self.gamma
                })

        nodes: List[Dict[List[bool], TimeslotState]] = []

        # nodes.append({
        #     tuple([modulate(False)]*self.reg_size): {
        #         'alpha': 0,
        #         'beta': 0,
        #         'gamma_to_next': 0,
        #         'transition_bit': None,
        #         'previous': None,
        #         'output': None
        #     }
        # })

        

        for _ in range(len(systematic_sequence)+1):
            timeslot = {}

            for state in product([False, True], repeat=self.reg_size):
                timeslot[tuple(state)] = TimeslotState(0, 0, {
                    False: 0,
                    True: 0
                })
            nodes.append(timeslot)

        nodes[0][tuple([False]*self.reg_size)].alpha = 1.0

        # The guarantee of all-zero final state
        if self.terminate == True:
            nodes[-1][tuple([False]*self.reg_size)].beta = 1.0
        else:
            for memory, node in nodes[-1].items():
                node.beta = 1/(2**self.reg_size)
        # Gamma computation
        t = 0
        for timeslot in nodes[:-1]:
            for memory, node in timeslot.items():
                for next_bit in [False, True]:
                    local_gamma = compute_gamma(systematic_sequence[t], parity_sequence[t], modulate(self.trellis[memory][next_bit]['output']))
                    
                    node.gamma[next_bit] = local_gamma
            t += 1

        # Forward recursion
        t = 1
        for timeslot in nodes[1:]:
            for memory, node in timeslot.items():
                previous_registers = self.trellis[memory]['previous_state']
                bit_to_state = self.trellis[memory]['bit_to_state']
                previous_node_1 = nodes[t-1][previous_registers[0]]
                previous_node_2 = nodes[t-1][previous_registers[1]]
                
                alpha = max_star(previous_node_1.alpha+previous_node_1.gamma[bit_to_state],  previous_node_2.alpha+previous_node_2.gamma[bit_to_state])

                node.alpha = alpha

            t += 1

        # Backward recursion
        t = len(nodes)-2
        for timeslot in nodes[-2::-1]:
            # print(t)
            for memory, node in timeslot.items():
                beta = 0

                for bit in [False, True]:
                    next_register = self.trellis[memory][bit]['next_state']
                    beta = max_star(nodes[t+1][next_register].beta+node.gamma[bit], beta)
                
                node.beta = beta

            t -= 1
        # Computation of Log Likelihood Ratios
        llrs = []
        t = 0
        for timeslot in nodes[:-1]:
            # max_1 = max_star()
            top = 0 
            bottom = 0
            for memory, node in timeslot.items():
                alpha = node.alpha
                beta_one = nodes[t+1][self.trellis[memory][True]['next_state']].beta
                beta_zero = nodes[t+1][self.trellis[memory][False]['next_state']].beta
                top += alpha * beta_one * node.gamma[True]
                bottom += alpha * beta_zero * node.gamma[False]

            llrs.append(np.log(top/bottom))

            t += 1

        # print(llrs)

        # print(nodes)

        decoded_sequence = [True if llr > 0 else False for llr in llrs]

        return decoded_sequence, llrs


def gates(register: List[bool]):
    return [register[0] ^ register[1] ^ register[2]]

import json
import time

if __name__ == '__main__':
    simulation_data = {}
    start = 0
    end = 0

    bers = []

    if sys.argv[1] == 'viterbi':
        rsc = RSC(gates, [1], 3, 2)

        BITS = 10000 * 100

        snr_in_db_range = np.arange(0, 9, 0.5)

        start = time.time()

        for index in range(len(snr_in_db_range)):
            local_bers = []

            snr_in_db = snr_in_db_range[index]

            simulation_data[snr_in_db] = {}

            snr_linear = 10**(snr_in_db/10)

            noise_spectral_dens = 1/snr_linear

            for i in range(10):
                simulation_data[snr_in_db][i] = {}
                
                uncoded = np.random.randint(0, 2, BITS, bool)

                msg = rsc.encode(uncoded)

                signal = 1 - 2*np.asarray(msg)

                print('n0', noise_spectral_dens)

                randoms = [box_muller() for _ in range(len(msg))]

                noise = sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                received = signal + noise

                demodulated = [1 if r < 0 else 0 for r in received]
                
                decoded = []

                
                decoded = rsc.decode(demodulated)

                n_errors = sum(decoded != uncoded)
                ber = n_errors/BITS

                simulation_data[snr_in_db][i]['number_of_errors'] = int(n_errors)
                simulation_data[snr_in_db][i]['bit_error_rate'] = float(ber)

                local_bers.append(ber)

            bers.append(average(local_bers))
        end = time.time()
        print('Elapsed: ', end - start)

        print(bers)

    elif sys.argv[1] == 'bcjr':
        rsc = RSC(gates, [1], 3, 2)

        BITS = int(sys.argv[2]) # 10000

        snr_in_db_range = np.arange(0, 9, 0.5)

        start = time.time()

        for index in range(len(snr_in_db_range)):
            local_bers = []

            snr_in_db = snr_in_db_range[index]

            simulation_data[snr_in_db] = {}

            snr_linear = 10**(snr_in_db/10)

            noise_spectral_dens = 1/snr_linear

            for i in range(int(sys.argv[3])): # 10):
                simulation_data[snr_in_db][i] = {}
                
                uncoded = np.random.randint(0, 2, BITS, bool)

                msg = rsc.encode(uncoded)

                signal = 1 - 2*np.asarray(msg)

                print('n0', noise_spectral_dens)

                randoms = [box_muller() for _ in range(len(msg))]

                noise = sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                received = signal + noise
                                
                decoded = rsc.bcjr_decode(received[0::2], received[1::2], variance=noise_spectral_dens/2)[0]

                n_errors = sum(decoded != uncoded)
                ber = n_errors/BITS

                simulation_data[snr_in_db][i]['number_of_errors'] = int(n_errors)
                simulation_data[snr_in_db][i]['bit_error_rate'] = float(ber)

                local_bers.append(ber)

            bers.append(average(local_bers))

        end = time.time()
        print('Elapsed: ', end - start)

        print(bers)

    

    elif sys.argv[1] == 'trellis':
        rsc = RSC(gates, [1], 3, 2, terminate=False)

        while True:

            random_bits = np.random.randint(0,2,2000,bool)
            # print(random_bits)



            coded = rsc.encode(random_bits)
            # print(coded)

            modulated_bits = modulate(coded)

            # print('Modulated: ', modulated_bits)
            

            systematic = modulated_bits[::2]
            parity = modulated_bits[1::2]

            output = rsc.bcjr_decode(systematic, parity)
            decoded = output[0]
            # print(output[1])

            # print(decoded)

            def pairwise_compare(a: list, b: list):
                for x, y in zip(a, b):
                    if x != y: return False

                return True

            if not pairwise_compare(random_bits, decoded):
                print('We don\'t have correct sequence!')

    if sys.argv[1] == 'bcjr' or sys.argv[1] == 'viterbi':
        if not os.path.exists('rsc'):
                os.makedirs('rsc')

        with open(f'rsc/rsc_{sys.argv[1]}_aggregated_bers_bits_{sys.argv[2]}_n_tests_{sys.argv[3]}_{end}.csv', mode='w') as file:
                writer = csv.writer(file)
                writer.writerows(map(lambda x: [x], bers))

        with open(f'rsc/rsc_{sys.argv[1]}_raw_data_bits_{sys.argv[2]}_n_tests_{sys.argv[3]}_{end}.json', mode='w') as file:
            json.dump(simulation_data, file)