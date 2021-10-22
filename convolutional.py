import os
import sys
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
import time
import json

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

    def _push_reg(self, bit: bool):
        last_bit = self.register[0]

        self.register = self.register[1:] + [bit]

        output = self.gates([last_bit] + self.register)

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

    def fill_trellis(self):
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
            self.trellis[previous['state']]['previous_state'].append(previous['previous_state'])



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

        # with open('./nodes.txt', 'w') as file:
        #     file.write(str(nodes))

        #   Backtracking the path with the lowest weight (NOTE: the survivor path)

        survivor_end = {
            'weight': Infinity
        }


        for state, node in nodes[-1].items():
            if node['weight'] < survivor_end['weight']:
                survivor_end = node

        # print(nodes[-1])

        # print(survivor_end)

        reversed_decoded: List[bool] = []

        while survivor_end.get('previous') is not None:
            reversed_decoded.append(survivor_end['transition_bit'])
            nodes.pop()
            survivor_end = nodes[-1][survivor_end['previous']]

        return reversed_decoded[::-1]

        
def gates_1(input: List[bool]):
    return [input[-1] ^ input[-3] ^ input[-4], input[-1] ^ input[-2] ^ input[-3]]

def gates_2(input: List[bool]):
    return [input[-1] ^ input[-2] ^ input[-3] ^ input[-4], input[-1] ^ input[-2] ^ input[-4]]

if __name__ == '__main__':
    if sys.argv[-1] == 'simulation':
        simulation_data = {}

        gates = []
        
        if sys.argv[-2] == 'gates_1':
            conv_code = ConvCode(gates_1, 3)
        elif sys.argv[-2] == 'gates_2':
            conv_code = ConvCode(gates_2, 3)

        BITS = int(sys.argv[1])

        snr_in_db_range = np.arange(0, 8, 0.5)

        bers = []

        start = time.time()

        for index in range(len(snr_in_db_range)):
            local_bers = []

            snr_in_db = snr_in_db_range[index]

            simulation_data[snr_in_db] = {}

            snr_linear = 10**(snr_in_db/10)

            # noise_spectral_dens = 1/snr_linear
            noise_variance = 0.5/snr_linear

            print('n0', noise_variance)

            for i in range(int(sys.argv[2])):
                simulation_data[snr_in_db][i] = {}
                
                uncoded = np.random.randint(0, 2, BITS, bool)

                msg = conv_code.encode(uncoded)

                signal = 1 - 2*np.asarray(msg)

                randoms = [box_muller() for _ in range(len(msg))]

                noise = sqrt(noise_variance)*np.asarray(randoms)

                received = signal + noise

                demodulated = [1 if r < 0 else 0 for r in received]
                
                decoded = []                
                decoded = conv_code.decode(demodulated)

                n_errors = sum(decoded != uncoded)
                ber = n_errors/BITS

                simulation_data[snr_in_db][i]['number_of_errors'] = int(n_errors)
                simulation_data[snr_in_db][i]['bit_error_rate'] = float(ber)

                local_bers.append(ber)


            bers.append(average(local_bers))
        end = time.time()
        print(bers)
        print('Elapsed: ', end - start)
        
        if not os.path.exists('convolutional'):
            os.makedirs('convolutional')

        with open(f'convolutional/convolutional_{sys.argv[-2]}_n_bits_{sys.argv[1]}_n_tests_{sys.argv[2]}_aggregated_bers_{time.strftime("%m-%d-%Y-%H-%M-%S", time.localtime(end))}.csv', mode='w') as file:
                writer = csv.writer(file)
                writer.writerows(map(lambda x: [x], bers))

        with open(f'convolutional/convolutional_{sys.argv[-2]}_n_bits_{sys.argv[1]}_n_tests_{sys.argv[2]}_raw_data_{time.strftime("%m-%d-%Y-%H-%M-%S", time.localtime(end))}.json', mode='w') as file:
            json.dump(simulation_data, file)
    
    elif sys.argv[-1] == 'trellis':
        code = ConvCode(gates_1, 3, 2)

        print(code.trellis)
