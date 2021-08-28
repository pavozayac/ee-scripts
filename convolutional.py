from typing import Callable, List
from itertools import product
import numpy as np

def chunks(lst: list, n: int):
       for i in range(0, len(lst), n):
              yield lst[i:i + n]


def push_shift(next: bool, arr: List[bool]):
    arr.pop[0]
    arr.append(next)

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
        output = self.gates(self.register)
        self.register = self.register[1:] + [bit]

        return output

    def fill_trellis(self):
        posregs = product([False, True], repeat=self.reg_size)

        for reg in posregs:
            self.trellis[reg] = {}
            self.register = list(reg).copy()

            for next_bit in [False, True]:
                output = self.push_reg(next_bit)

                self.trellis[reg][next_bit] = {
                    'next_state': tuple(self.register),
                    'output': output
                }


    def encode(self, bits: List[int]) -> List[int]:
        boolean_bits = map(lambda x: bool(x), bits)
        output: List[int] = []

        for bit in boolean_bits:
            output += map(lambda x: int(x), self.push_reg(bit))

        return output

    def decode(self, bits: List[int]):
        boolean_bits = list(map(lambda x: bool(x), bits))
        grouped_bits = chunks(boolean_bits, self.out_size)







def half_gates(input: List[bool]):
    return [input[0] ^ input[1] ^ input[2], input[1]]

if __name__ == '__main__':
    b = ConvCode(half_gates, reg_size=3, out_size=2)

    test = np.random.randint(0, 2, 10000 * 100, int)

    print(len(b.encode(test)))