from typing import List
import random
import numpy as np
import os
import csv

class Interleaver:
    def __init__(self, length: int):
        self.size = length

        if os.path.isfile(f'interleaver_{length}.csv'):
            self.indices = []
            with open(f'interleaver_{length}.csv', newline='\n') as file:
                reader = csv.reader(file, delimiter=' ')

                for row in reader:
                    # print(row)
                    self.indices.append(int(row[0]))

        else:
            self.size = length
            s = list(range(0, length))
            self.indices = []

            for _ in range(length):
                choice = random.choice(s)
                s.remove(choice)
                self.indices.append(choice)

            with open(f'interleaver_{length}.csv', mode='w') as file:
                writer = csv.writer(file)

                writer.writerows(map(lambda x: [x], self.indices))

    def interleave(self, sequence: list) -> list:
        output = []

        for index in self.indices:
            output.append(sequence[index])

        return output

    def deinterleave(self, sequence: list) -> list:
        output = [None]*self.size

        for number, index in enumerate(self.indices):
            output[index] = sequence[number]

        return output
            




if __name__ == '__main__':
    i = Interleaver(10)
    print(i.indices)

    random_bits = np.random.randint(0, 10, 10, int)
    print(random_bits)

    interleaved = i.interleave(random_bits)
    print(interleaved)

    deinterleaved = i.deinterleave(interleaved)
    print(deinterleaved)
    

