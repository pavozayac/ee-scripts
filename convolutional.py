import os
from typing import Callable, Dict, List
from itertools import product
import numpy as np
from numpy.core.numeric import Infinity
from numpy.lib.function_base import append
from numpy.lib.function_base import average
import time
from box_muller import box_muller
from math import sqrt
import csv

def chunks(lst: list, n: int):
       for i in range(0, len(lst), n):
            yield lst[i:i + n]


def push_shift(next: bool, arr: List[bool]):
    arr.pop[0]
    arr.append(next)

def hamming_distace(a: List[bool], b: List[bool]):
    
    c = [x ^ y for x, y in zip(a, b)]

    return sum(c)
    

class ConvCode:
    trellis: dict
    register_size: int
    output_size: int
    register: List[bool]
    gates: Callable[[List[bool]], List[bool]]

    
    def __init__(self, gates, reg_size = 3, out_size = 2, ):
        self.trellis = {}
        self.reg_size = reg_size
        self.out_size = out_size
        self.gates = gates
        self.register = [False]*reg_size

        self.fill_trellis()

    def push_reg(self, bit: bool):
        self.register = self.register[1:] + [bit]

        output = self.gates(self.register)


        return output

    def fill_trellis(self):
        posregs = product([False, True], repeat=self.reg_size)

        for reg in posregs:
            self.trellis[reg] = {}

            for next_bit in [False, True]:
                self.register = list(reg).copy()

                output = self.push_reg(next_bit)

                self.trellis[reg][next_bit] = {
                    'next_state': tuple(self.register),
                    'output': output
                }


    def encode(self, bits: List[bool]) -> List[bool]:
        self.register = [False]*self.reg_size

        output: List[int] = []

        for bit in bits:
            output += self.push_reg(bit)


        return output

    def decode(self, bits: List[bool]):
        grouped_bits = chunks(bits, self.out_size)

        #   The nodes is a list of time arranged lists of all possible states of the register with assigned cumulative weights
        nodes = []

        nodes.append({
            tuple([False]*self.reg_size): {
                'weight': 0,
                'transition_bit': None,
                'previous': None,
                'output': None
            }
        })

        #   Generating the weights of the nodes
        for received_output in grouped_bits:
            next_timeslot = {}

            for state, node in nodes[-1].items():
                for next_bit in [False, True]:
                    next_trellis = self.trellis[state][next_bit]
                    distance = hamming_distace(next_trellis['output'], received_output)

                    if next_timeslot.get(next_trellis['next_state']) is not None:
                        if distance + node['weight'] < next_timeslot[next_trellis['next_state']]['weight']:
                            next_timeslot[next_trellis['next_state']] = {
                                'weight': distance + node['weight'],
                                'transition_bit': next_bit,
                                'previous': state,
                                'output': next_trellis['output']
                            }
                    else:
                        next_timeslot[next_trellis['next_state']] = {
                            'weight': distance + node['weight'],
                            'transition_bit': next_bit,
                            'previous': state,
                            'output': next_trellis['output']
                        }

            nodes.append(next_timeslot)

        with open('./nodes.txt', 'w') as file:
            file.write(str(nodes))

        #   Backtracking the path with the lowest weight (NOTE: the survivor path)

        survivor_end = {
            'weight': Infinity
        }


        for state, node in nodes[-1].items():
            if node['weight'] < survivor_end['weight']:
                survivor_end = node

        print(nodes[-1])

        print(survivor_end)

        reversed_decoded: List[bool] = []

        while survivor_end.get('previous') is not None:
            reversed_decoded.append(survivor_end['transition_bit'])
            nodes.pop()
            survivor_end = nodes[-1][survivor_end['previous']]

        return reversed_decoded[::-1]

        


def half_gates(input: List[bool]):
    return [input[0] ^ input[1] ^ input[2], input[2]]

if __name__ == '__main__':

    gates = []
        
    conv_code = ConvCode(half_gates)

    BITS = 10000 * 100

    snr_in_db_range = np.arange(0, 9, 0.5)

    bers = []

    start = time.time()

    for index in range(len(snr_in_db_range)):
            local_bers = []

            snr_in_db = snr_in_db_range[index]

            snr_linear = 10**(snr_in_db/10)

            noise_spectral_dens = 1/snr_linear

            for i in range(10):
                    
                    
                    uncoded = np.random.randint(0, 2, BITS, bool)

                    msg = conv_code.encode(uncoded)

                    signal = 1 - 2*np.asarray(msg)

                    print('n0', noise_spectral_dens)

                    randoms = [box_muller() for _ in range(len(msg))]

                    noise = sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                    received = signal + noise

                    demodulated = [1 if r < 0 else 0 for r in received]
                    
                    decoded = []

                    
                    decoded = conv_code.decode(demodulated)

                    n_errors = sum(decoded != uncoded)
                    ber = n_errors/BITS
                    local_bers.append(ber)


            bers.append(average(local_bers))
    end = time.time()
    print('Elapsed: ', end - start)
    
    with open('convolutional.csv', mode='w') as file:
            writer = csv.writer(file)
            writer.writerows(map(lambda x: [x], bers))
