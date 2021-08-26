from typing import List


def hamming_encode (arr: list):
    if len(arr) % 4 != 0:
        raise Exception("Message not aligned to 4 bits")
    
    blocks = []
    for i in range(int(len(arr)/4)):
        blocks.append(arr[i*4:i*4+4])

    for n in blocks:
        m = list(n)
        m.insert(0, 0)
        m.insert(1, 0)
        m.insert(3, 0)
        m[0] = (m[0]+m[2]+m[4]+m[6])%2
        m[1] = (m[1]+m[2]+m[5]+m[6])%2
        m[3] = (sum(m[3:7]))%2

    return [bit for block in blocks for bit in block]

def hamming_correct(arr: list):
    blocks = []
    for i in range(int(len(arr)/7)):
        blocks.append(arr[i*7:i*7+7])

    for n in blocks:
        m = list(n)
        syndrome = []
        syndrome.append((sum(m[3:7]))%2)
        syndrome.append((m[1]+m[2]+m[5]+m[6])%2)
        syndrome.append((m[0]+m[2]+m[4]+m[6])%2)
        if syndrome != [0,0,0]:
            position = -1
            for bit in syndrome:
                position = (position<<1) | bit
            if m[position] == 0:
                m[position] = 1
            else:
                m[position] = 0
            
    return [bit for block in blocks for bit in block]

def hamming_decode(arr: list):
    
    blocks = []
    for i in range(int(len(arr)/7)):
        blocks.append(arr[i*7:i*7+7])

    output = []
    for m in blocks:
        output.append([m[2]]+m[4:7])

    return [bit for block in output for bit in block]

    
if __name__ == '__main__':
    enc = hamming_encode([1,0,1,1,1,0,0,1,0,0,1,0])
    enc[0] = 1
    print(enc)
    correct = hamming_correct(enc)
    print(correct)
    print(hamming_decode(correct))