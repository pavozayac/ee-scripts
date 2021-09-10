import random
from typing import List
import numpy as np

class Interleaver:
    def __init__(self, length: int):
        self.size = length
        s = list(range(0, length))
        self.indices = []

        for _ in range(length):
            choice = random.choice(s)
            s.remove(choice)
            self.indices.append(choice)

    def interleave(self, sequence: List[bool]) -> List[bool]:
        output = []

        for index in self.indices:
            output.append(sequence[index])

        return output

    def deinterleave(self, sequence: List[bool]) -> List[bool]:
        output = [None]*self.size

        for number, index in enumerate(self.indices):
            output[index] = sequence[number]

        return output
            




if __name__ == '__main__':
    i = Interleaver(10)
    print(i.indices)

    random = np.random.randint(0, 10, 10, int)
    print(random)

    interleaved = i.interleave(random)
    print(interleaved)
    deinterleaved = i.deinterleave(interleaved)
    print(deinterleaved)
    

