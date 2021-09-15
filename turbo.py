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
from rsc import RSC

class TurboCode():
    code1: RSC
    code2: RSC
    interleaver: Interleaver

    def __init__(self, code1, code2, interleaver):
        self.code1 = code1
        self.code2 = code2
        self.interleaver = interleaver

    def encode(self, sequence: List[bool]) -> List[bool]:
        output1 = chunks(self.code1.encode(sequence), 2)

        interleaved = self.interleaver.interleave(sequence)
        output2 = self.code2.encode(interleaved)[1::2]
        print('O2: ', output2)
        final = list(zip(output1, output2))

        return final

    def decode(self, sequence: List[float], max_iterations = 5):
        systematic = sequence[0::3]
        parity1 = sequence[1::3]
        parity2 = sequence[2::3]
        parity2_deinterleaved = self.interleaver.deinterleave(parity2)



if __name__ == '__main__':
    # random_bits = np.random.randint(0, 2, 10, bool)
    # print(random_bits)

    # def gates(input: List[bool]):
    #     return [input[0] ^ input[2]]

    # rsc = RSC(gates, [1], 3, 2)

    # coded = rsc.encode(random_bits)
    # decoded = rsc.decode(coded)
    # print(coded)
    # print(decoded)

    def gates(register: List[bool]):
        return [register[0] ^ register[1] ^ register[2]]

    rsc1 = RSC(gates, [1], 3, 2)
    rsc2 = RSC(gates, [1], 3, 2)
    interleaver = Interleaver(128)

    code = TurboCode(rsc1, rsc2, interleaver)
    random_bits = np.random.randint(0,2,10,bool)

    coded = rsc1.encode(random_bits)
    print(len(coded))
    decoded = rsc1.decode(random_bits)

    print(coded)
    print(decoded)