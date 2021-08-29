import os
from typing import Callable, Dict, List
from itertools import product
import numpy as np
from numpy.core.numeric import Infinity
from numpy.lib.function_base import append

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


    def encode(self, bits: List[int]) -> List[int]:
        self.register = [False]*self.reg_size

        boolean_bits = map(lambda x: bool(x), bits)
        output: List[int] = []

        for bit in boolean_bits:
            output += map(lambda x: int(x), self.push_reg(bit))


        return output

    def decode(self, bits: List[int]):
        boolean_bits = list(map(lambda x: bool(x), bits))
        grouped_bits = chunks(boolean_bits, self.out_size)

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
    b = ConvCode(half_gates, reg_size=3, out_size=2)

    print(b.trellis)

    test = np.random.randint(0, 2, 10000, int)

    encoded = b.encode(test)

    print(len(encoded))

    decoded = b.decode(encoded)

    print(len(decoded))

    def check_equivalence(a: list, b: list):
        ind = 0

     
        c = []

        for index, (x, y) in enumerate(zip(a, b)):
            if x != y:
                c.append(1)
                ind = index
            else:
                c.append(0)


        return sum(c), ind

    print(check_equivalence(test, decoded))


{(False, False, False): {False: {'next_state': (False, False, False), 'output': [False, False]}, True: {'next_state': (False, False, True), 'output': [False, False]}}, (False, False, True): {False: {'next_state': (False, True, False), 'output': [True, True]}, True: {'next_state': (False, True, True), 'output': [True, True]}}, (False, True, False): {False: {'next_state': (True, False, False), 'output': [True, False]}, True: {'next_state': (True, False, True), 'output': [True, False]}}, (False, True, True): {False: {'next_state': (True, True, False), 'output': [False, True]}, True: {'next_state': (True, True, True), 'output': [False, True]}}, (True, False, False): {False: {'next_state': (False, False, False), 'output': [True, False]}, True: {'next_state': (False, False, True), 'output': 
[True, False]}}, (True, False, True): {False: {'next_state': (False, True, False), 'output': [False, True]}, True: {'next_state': (False, True, True), 'output': [False, True]}}, (True, True, False): 
{False: {'next_state': (True, False, False), 'output': [False, False]}, True: {'next_state': (True, False, True), 'output': [False, False]}}, (True, True, True): {False: {'next_state': (True, True, False), 'output': [True, True]}, True: {'next_state': (True, True, True), 'output': [True, True]}}}