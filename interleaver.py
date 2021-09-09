import random
from typing import List

class Interleaver:
    def __init__(self, length: int):
        self.s = list(range(0, length))
        self.indices = []

        for i in range(length):
            choice = random.choice(self.s)
            self.s.remove(choice)
            self.indices.append(choice)

    def interleave(sequence: List[bool]):
        output = []

        for bit in sequence:
            




if __name__ == '__main__':
    i = Interleaver(10)
    print(i.indices)
    

