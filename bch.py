import numpy as np 
import galois

class BCH():
    
    def __init__(self, order: int, irreducible: galois.Poly, generator: list, t) -> None:
        self.field = galois.GF(2**order, irreducible)
        self.order = order
        self.field.display('power')
        self.generator = galois.Poly(generator, self.field)
        #t is the error-correcting capacity of the code
        self.t = t

    def encode(self, data: list) -> list:

        msg = galois.Poly(data, self.field)
        codeword: galois.Poly = self.generator * msg

        encoded = [int(b) for b in codeword.coeffs]
        padded = [0]*(127-len(encoded)) + encoded
        return padded

    def correct(self, data: list) -> list:
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
                print(loc_coeffs.flatten())
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

    def decode(self, data: list) -> list:
        datapoly = galois.Poly(data, self.field)
        decoded_poly = datapoly/self.generator

        decoded_data = [int(b) for b in decoded_poly.coeffs]
        padded = [0]*(113-len(decoded_data)) + decoded_data
        return padded
    

if __name__ == '__main__':
    bch127 = BCH(7, galois.Poly([1, 0, 0, 0, 1, 0, 0, 1], galois.GF(2)), [1,0,0,0,0,1,1,0,1,1,1,0,1,1,1], 2)

    msg = np.random.randint(2, size=113)
    codeword = bch127.encode(msg)
    modified = codeword.copy()
    print(codeword)
    modified[17] = 1 if modified[17] == 0 else 0
    #codeword[30] = 1 if codeword[30] == 0 else 0

    corrected = bch127.correct(modified)
    #print(modified)
    #print(len(bch127.correct(modified)))
    decoded = bch127.decode(corrected)
    print(len(decoded))
    #decoded[0] = 1 if decoded[0] == 0 else 0
    print(decoded == msg)