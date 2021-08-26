import numpy as np 
import galois
from functools import lru_cache

class BCH():
    
    def __init__(self, order: int, irreducible: galois.Poly, generator: list, t) -> None:
        self.field = galois.GF(2**order, irreducible)
        self.order = order
        self.field.display('power')
        self.generator = galois.Poly(generator, self.field)
        #t is the error-correcting capacity of the code
        self.t = t

    @lru_cache(maxsize=None)
    def encode(self, data: tuple) -> list:
        data = list(data)

        msg = galois.Poly(data, self.field)
        codeword: galois.Poly = self.generator * msg

        encoded = [int(b) for b in codeword.coeffs]
        padded = [0]*(15-len(encoded)) + encoded
        return padded
    @lru_cache(maxsize=None)
    def correct(self, data: tuple) -> list:
        data = list(data)

        received = galois.Poly(data, self.field)

        syndromes = []

        for i in range(2*self.t):
            syndromes.append(received(self.generator.roots()[0]**(i+1)))

        for v in range(self.t, 0, -1):
            m = []

            #galois.log.log_naive

            for i in range(v):
                m.append(syndromes[i:v+i])

            if np.linalg.det(self.field(m)) != 0:
                svs = self.field([[x] for x in syndromes[v:2*v]])

                loc_coeffs = np.linalg.inv(self.field(m)) @ self.field(svs)
                #print(loc_coeffs.flatten())
                positions = galois.Poly([*loc_coeffs.flatten(), 1], self.field).roots()
                for p in positions:
                    i = -1
                    temp = p
                    while temp != 1:
                        temp = temp/self.field.primitive_element
                        i += 1
                    
                    data[i] = 1 if data[i] == 0 else 0

                #print(self.field.)

                #for loc in p:
                #   print(loc[0])
                return data
            return data

    @lru_cache(maxsize=None)
    def decode(self, data: tuple) -> list:
        data = list(data)

        datapoly = galois.Poly(data, self.field)
        decoded_poly = datapoly/self.generator

        decoded_data = [int(b) for b in decoded_poly.coeffs]
        padded = [0]*(7-len(decoded_data)) + decoded_data
        return padded

def bch15_7():
    return BCH(4, galois.Poly([1, 0, 0, 1, 1], galois.GF(2)), [1,1,1,0,1,0,0,0,1], 2)

if __name__ == '__main__':
    bch15_7 = BCH(4, galois.Poly([1, 0, 0, 1, 1], galois.GF(2)), [1,1,1,0,1,0,0,0,1], 2)

    true = 0
    total = 0

    for i in range(0,15):
        for j in range(0,15):

            if i == j:
                continue

            total += 1


            msg = np.random.randint(2, size=7)
            codeword = bch15_7.encode(msg)
            modified = codeword.copy()
            print(codeword)
            modified[i] = 1 if modified[i] == 0 else 0
            modified[j] = 1 if modified[j] == 0 else 0
            #codeword[30] = 1 if codeword[30] == 0 else 0

            corrected = bch15_7.correct(modified)
            #print(modified)
            #print(len(bch127.correct(modified)))
            decoded = bch15_7.decode(corrected)
            #print(len(decoded))
            #decoded[0] = 1 if decoded[0] == 0 else 0
            #print(decoded == msg)
            if all(decoded == msg):
                true += 1

    print(true, total)