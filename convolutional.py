from typing import Callable, List

def push_shift(num: int, arr: list):
    for i in range(len(arr)-1, 0, -1):
        arr[i] = arr[i-1]
    arr[0] = num

class ConvCode:
    trellis: dict
    register_size: int
    output_size: int
    register: list
    gates: Callable[[List[int]], List[int]]

    
    def __init__(self, gates, reg_size = 3, out_size = 2, ):
        self.trellis = {}
        self.reg_size = reg_size
        self.out_size = out_size
        self.gates = gates
        self.register = [0]*reg_size

    def push_reg(self, bit: int):
        output = self.gates(self.register)
        self.register.pop(0)
        self.register.append(bit)

        return output

    def generate_trellis(self):

        for i in range(2):
            for j in range(2):
                

                self.trellis[[i, j]] = {

                }

    def encode(self, bits) -> list:

        encoded_seq = []

        registers = [0, 0]

        for bit in bits:
            o1 = (bit + registers[1]) % 2
            o2 = (bit + registers[0] + registers[1]) % 2

            encoded_seq.append(o1)
            encoded_seq.append(o2)

            push_shift(bit, registers)

        return encoded_seq


if __name__ == '__main__':
    
    b = ConvCode()

    msg = [1,1,0]

    print(b.encode(msg))