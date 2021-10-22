import math
from numpy import random
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

def flatten(list_like):
    if isinstance(list_like, (list, tuple, set, range)):
        for sub in list_like:
            yield from flatten(sub)
    else:
        yield list_like

def pairwise_compare(a: list, b: list):
    for x, y in zip(a, b):
        if x != y: return False

    return True

class RSCMode(Enum):
    SYSTEMATIC = 'systematic'
    NONSYSTEMATIC = 'nonsystematic'

class RSCDecoder(Enum):
    BCJR = 'bcjr'
    VITERBI = 'viterbi'

def modulate(a: bool) -> float:
    if isinstance(a, bool):
        return 2*float(a)-1
    else:
        return [2*float(b)-1 for b in a]


def max_star(values: list):
    values = values.copy()
    buffer = values[0]
    for value in values[1:]:
        if value == float('-inf'):
            buffer = buffer
        else:
            buffer = max(buffer, value) + np.log(1 + math.e**(-np.abs(value-buffer)))

    return buffer
    


class RSC(ConvCode):
    # RSC should always have one gate less than its output

    def __init__(self, gates, recursive_indices: List[int], reg_size = 3, out_size = 2, terminate = False, decoder = RSCDecoder.VITERBI):
        self.trellis = {}
        self.reg_size = reg_size
        self.out_size = out_size
        self.gates = gates
        self.register = [False]*reg_size
        self.recursive_indices = recursive_indices
        self.terminate = terminate
        self.decoder = decoder
        self.snr = None

        #if decoder == RSCDecoder.VITERBI:
        self.fill_trellis_a()
        #else:
            #self.fill_trellis_bcjr()

    def set_snr(self, snr):
        self.snr = snr

    def fill_trellis_bcjr(self):
        posregs = product([modulate(False), modulate(True)], repeat=self.reg_size)

        for reg in posregs:
            prepared = reg + [False]*3
            


    def fill_trellis_a(self):
        posregs = product([False, True], repeat=self.reg_size)

        previous_to_fill = []

        for reg in posregs:
            self.trellis[reg] = {}
            self.trellis[reg]['previous_state'] = []
            self.trellis[reg]['bit_to_state'] = reg[-1]

            for next_bit in [False, True]:
                self.register = list(reg).copy()

                output = self._push_reg(next_bit)

                self.trellis[reg][next_bit] = {
                    'next_state': tuple(self.register),
                    'output': output
                }

                previous_to_fill.append({
                    'state': tuple(self.register),
                    'previous_bit': next_bit,
                    'previous_state': reg
                })

        for previous in previous_to_fill:
            self.trellis[previous['state']]['previous_state'].append({
                    'register': previous['previous_state'],
                    'bit_to_state': previous['previous_bit']
                }
            )

        print(self.trellis)
    
    def _push_reg(self, bit: bool):
        recursive_bit = 0

        last_bit = self.register[0]

        for i in self.recursive_indices:
            # NOTE The first input is always added to the recursive registers, so the argument in RSC() cannot be the full feedback polynomial!!!
            recursive_bit = recursive_bit ^ ([last_bit]+self.register+[bit])[i]

        self.register = self.register[1:] + [recursive_bit]

        output = [bit] + self.gates([last_bit] + self.register)

        return output

    def push_reg(self, bit: bool):
        # recursive_bit = bit
        
        # for i in self.recursive_indices:
        #     recursive_bit = recursive_bit ^ self.register[i]

        # self.register = self.register[1:] + [recursive_bit]
        # output = [bit] + self.gates(self.register)
        output = self.trellis[tuple(self.register)][bit]['output']
        self.register = self.trellis[tuple(self.register)][bit]['next_state']
        return output


    def encode(self, bits: List[bool]) -> List[bool]:
        self.register = [False]*self.reg_size

        output: List[bool] = []

        for bit in bits:
            output += self.push_reg(bit)

        if self.terminate == True:
            for _ in range(self.reg_size):
                termination_bit = self.register[self.recursive_indices[0]]

                for index in self.recursive_indices[1:]:
                    termination_bit = termination_bit ^ self.register[index]

                output += self.push_reg(termination_bit)
                print(termination_bit)
                print(self.register)


        # print(self.register)

        return output

    #This gamma is defined specifically for BPSK modulation
    def compute_gamma(self, signal: float, parity: float, prototype_output: List[float], a_priori_log = None) -> float:
        if a_priori_log:
            # print((prototype_output[1]*a_priori_log)/2)
            # print('we have a priori', (0.5)*prototype_output[0]*a_priori_log)
            return (self.Lc/2)*(prototype_output[0]*signal+prototype_output[1]*parity) + (prototype_output[0]*a_priori_log)/2
        else: 
            return (self.Lc/2)*(prototype_output[0]*signal+prototype_output[1]*parity)

    
    #takes a list of floating-point numbers and returns the log likelihood ratio, optionally takes a list of previous LLRs
    def bcjr_decode(self, systematic_sequence: List[float], parity_sequence: List[float], Le: Optional[List[float]] = None) -> List[float]:
        self.Lc = 2*self.snr

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

        for _ in range(len(systematic_sequence)+1):
            timeslot = {}

            for state in product([False, True], repeat=self.reg_size):
                # print(state)
                timeslot[tuple(state)] = TimeslotState(float('-inf'), float('-inf'), {
                    False: 0,
                    True: 0
                })
            nodes.append(timeslot)

        nodes[0][tuple([False]*self.reg_size)].alpha = np.log(1.0)

        # The guarantee of all-zero final state
        if self.terminate == True:
            nodes[-1][tuple([False]*self.reg_size)].beta = np.log(1.0)
        else:
            for memory, node in nodes[-1].items():
                node.beta = np.log(1/(2**self.reg_size))

        # Gamma computation
        t = 0
        for timeslot in nodes[:-1]:
            for memory, node in timeslot.items():
                for next_bit in [False, True]:
                    #print(self.trellis[memory][next_bit]['output'], modulate(self.trellis[memory][next_bit]['output']))
                    local_gamma = self.compute_gamma(
                        systematic_sequence[t], parity_sequence[t], 
                        modulate(self.trellis[memory][next_bit]['output']), 
                        a_priori_log=Le[t] if Le else None
                    )
                    
                    timeslot[memory].gamma[next_bit] = local_gamma
            t += 1

        # Forward recursion
        t = 1
        for timeslot in nodes[1:]:
            for memory, node in timeslot.items():
                previous_states = self.trellis[memory]['previous_state']
                #bit_to_state = self.trellis[memory]['bit_to_state']
                alphas = []
                for previous in previous_states:
                    prev_node = nodes[t-1][previous['register']]
                    bit_to_state = previous['bit_to_state']
                    
                    alphas.append(prev_node.alpha+prev_node.gamma[bit_to_state])
                    
                alpha = max_star(alphas)

                timeslot[memory].alpha = alpha

            t += 1

        # Backward recursion
        for t in range(len(nodes)-2,0,-1):
            timeslot = nodes[t]

            for memory in timeslot:
                betas = []

                for next_beta_bit in [False, True]:
                    next_register = self.trellis[memory][next_beta_bit]['next_state']
                    beta_component = nodes[t+1][next_register].beta
                    betas.append(beta_component+timeslot[memory].gamma[next_beta_bit])

                timeslot[memory].beta = max_star(betas)
        
        llrs = []
        for time in range(1,len(nodes)):
            # print(time)
            timeslot = nodes[time]

            ones = []
            zeros = []

            for memory, node in timeslot.items():                       
                beta_local = timeslot[memory].beta
                # to_state = self.trellis[memory]['bit_to_state'] 
                for previous in self.trellis[memory]['previous_state']:
                    alpha1 = nodes[time-1][previous['register']].alpha
                    gamma1 = nodes[time-1][previous['register']].gamma[previous['bit_to_state']]
                
                    if previous['bit_to_state'] == True:
                        ones.append(alpha1 + beta_local + gamma1)
                    else:
                        zeros.append(alpha1 + beta_local + gamma1)
            one = max_star(ones)
            zero = max_star(zeros)
            llrs.append(one-zero)
            
        for llr in llrs:
            if llr == 0: print('wth')

        decoded_sequence = [True if llr > 0 else False for llr in llrs]
     

        if Le == None:
            extrinsic_llrs = [(computed_llr - self.Lc*rec_sys) for computed_llr, rec_sys in zip(llrs, systematic_sequence)]
        else:
            extrinsic_llrs = [(computed_llr - self.Lc*rec_sys - ext) for computed_llr, rec_sys, ext in zip(llrs, systematic_sequence, Le)]
        # print('Le: ', extrinsic_llrs)
        # print('LLR: ', llrs)

        return decoded_sequence, extrinsic_llrs, llrs


def gates_1(register: List[bool]):
    return [register[-1] ^ register[-3] ^ register[-4]]

def gates_2(register: List[bool]):
    return [register[-1] ^ register[-2] ^ register[-3] ^ register[-4]]

import json
import time

if __name__ == '__main__':
    simulation_data = {}
    start = 0
    end = 0

    bers = []

    if sys.argv[-1] == 'gates_1':
        rsc = RSC(gates_1, [-1, -2, -3], 3, 2)
    elif sys.argv[-1] == 'gates_2':
        rsc = RSC(gates_2, [-1, -2, -4], 3, 2)


    if sys.argv[1] == 'viterbi':
        # NOTE !!! The first input is always a part of the recursive encoder, so it CANNOT be here!!!

        BITS = int(sys.argv[2])

        snr_in_db_range = np.arange(0, 8, 0.5)

        start = time.time()

        for index in range(len(snr_in_db_range)):
            local_bers = []

            snr_in_db = snr_in_db_range[index]

            simulation_data[snr_in_db] = {}

            snr_linear = 10**(snr_in_db/10)

            # noise_spectral_dens = 1/snr_linear
            noise_variance = 0.5/snr_linear

            print('n0', noise_variance)

            for i in range(int(sys.argv[3])):
                simulation_data[snr_in_db][i] = {}
                
                uncoded = np.random.randint(0, 2, BITS, bool)

                msg = rsc.encode(uncoded)

                signal = 1 - 2*np.asarray(msg)

                randoms = [box_muller() for _ in range(len(msg))]

                noise = sqrt(noise_variance)*np.asarray(randoms)

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
        BITS = int(sys.argv[2]) # 10000

        snr_in_db_range = np.arange(0, 8, 0.5)

        start = time.time()

        for index in range(len(snr_in_db_range)):
            local_bers = []

            snr_in_db = snr_in_db_range[index]

            simulation_data[snr_in_db] = {}

            snr_linear = 10**(snr_in_db/10)

            # noise_spectral_dens = 1/snr_linear
            noise_variance = 0.5/snr_linear

            print('n0', noise_variance)

            for i in range(int(sys.argv[3])): # 10):
                simulation_data[snr_in_db][i] = {}
                
                uncoded = np.random.randint(0, 2, BITS, bool)

                msg = rsc.encode(uncoded)

                signal = modulate(msg)


                randoms = [box_muller() for _ in range(len(msg))]

                noise = sqrt(noise_variance)*np.asarray(randoms)

                received = signal + noise

                rsc.set_snr(snr_linear)
                                
                decoded = rsc.bcjr_decode(received[0::2], received[1::2])[0]

                n_errors = sum(1 if a != b else 0 for a, b in zip(decoded, uncoded))
                ber = n_errors/BITS

                simulation_data[snr_in_db][i]['number_of_errors'] = int(n_errors)
                simulation_data[snr_in_db][i]['bit_error_rate'] = float(ber)

                local_bers.append(ber)

            bers.append(average(local_bers))

        end = time.time()
        print('Elapsed: ', end - start)

        print(bers)

    

    elif sys.argv[1] == 'trellis':
        rsc = RSC(gates_1, [1, 2], 3, 2, terminate=False)

        # print(rsc.trellis)
        #while True:

        random_bits = np.random.randint(0,2,10,bool)



        coded = rsc.encode(random_bits)
        print(coded, len(coded))
        print('Reg:', rsc.register)

        modulated_bits = modulate(coded)
        #modulated_bits[2] -= 2

        # print('Modulated: ', modulated_bits)
        

        systematic = modulated_bits[::2]
        parity = modulated_bits[1::2]

        output = rsc.bcjr_decode(systematic, parity)
        decoded = output[0]
        # print(output[1])

        # print(decoded)

        

        if not pairwise_compare(random_bits, decoded):
            print('We don\'t have correct sequence!')
        print(random_bits)
        print(coded)
        print(decoded)
        print(modulated_bits)


    if sys.argv[1] == 'bcjr' or sys.argv[1] == 'viterbi':
        if not os.path.exists('rsc'):
                os.makedirs('rsc')

        with open(f'rsc/rsc_{sys.argv[-1]}_{sys.argv[1]}_aggregated_bers_bits_{sys.argv[2]}_n_tests_{sys.argv[3]}_{time.strftime("%m-%d-%Y-%H-%M-%S", time.localtime(end))}.csv', mode='w') as file:
                writer = csv.writer(file)
                writer.writerows(map(lambda x: [x], bers))

        with open(f'rsc/rsc_{sys.argv[-1]}_{sys.argv[1]}_raw_data_bits_{sys.argv[2]}_n_tests_{sys.argv[3]}_{time.strftime("%m-%d-%Y-%H-%M-%S", time.localtime(end))}.json', mode='w') as file:
            json.dump(simulation_data, file)