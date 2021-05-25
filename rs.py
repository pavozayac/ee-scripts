from math import log2
from galois.field.poly_functions import primitive_element
import numpy as np 
import galois
from numpy.core.defchararray import index
from numpy.core.function_base import logspace

class RS():
    
    def __init__(self, n: int, k: int) -> None:
        self.n = n
        self.order = int(log2(n+1))
        self.t = int((n-k)/2)


        self.field = galois.GF(2**self.order, galois.conway_poly(2, self.order))
        self.field.display('power')

        a = self.field.primitive_element
        self.generator = galois.Poly.Roots(roots=[a**i for i in range(1, 2*self.t+1)], field=self.field)
        #print(self.generator)
        #t is the error-correcting capacity of the code

    def encode(self, data: list) -> list:

        symbols = []

        for i in range(int(len(data)/self.order)):
            symbols.append(data[i*self.order:i*self.order+self.order])

        coeffs = []

        for sym in symbols:
            coeffs.append(self.field.Elements()[int(''.join(str(bit) for bit in sym[::-1]), 2)])

        #print(coeffs)
        msg_poly = galois.Poly(coeffs[::-1], self.field)

        codeword = self.generator * msg_poly
        print('Codeword: ', codeword)

        binstr_coeffs = [bin(int(coeff))[2:] for coeff in codeword.coeffs]

        for index in range(len(binstr_coeffs)):
            binstr_coeffs[index] = [int(bit) for bit in binstr_coeffs[index]]
            binstr_coeffs[index] = (self.order-len(binstr_coeffs[index]))*[0] + binstr_coeffs[index]

        
        return [bit for sym in binstr_coeffs for bit in sym]

        '''
        msg = galois.Poly(data, self.field)
        codeword: galois.Poly = self.generator * msg

        encoded = [int(b) for b in codeword.coeffs]
        padded = [0]*(127-len(encoded)) + encoded
        return padded

        '''
    def correct(self, data: list) -> list:
        symbols = []

        for i in range(int(len(data)/self.order)):
            symbols.append(data[i*self.order:i*self.order+self.order])

        coeffs = []

        for sym in symbols:
            coeffs.append(self.field.Elements()[int(''.join(str(bit) for bit in sym), 2)])

        #print(coeffs)
        received_poly = galois.Poly(coeffs, self.field)
        print('Received: ', received_poly)

        #received = galois.Poly(data, self.field)

        syndromes = []

        for i in range(1, 2*self.t+1):
            syndromes.append(received_poly(self.field.primitive_element**i))

        print('syndromes: ', syndromes)
        for v in range(self.t, 0, -1):
            m = []

            #galois.log.log_naive

            for i in range(v):
                m.append(syndromes[i:v+i])

            if np.linalg.det(self.field(m)) != 0:
                svs = self.field([[x] for x in syndromes[v:2*v]])

                loc_coeffs = np.linalg.inv(self.field(m)) @ self.field(svs)
                #print(loc_coeffs.flatten())
                locator_poly = galois.Poly([*loc_coeffs.flatten(), 1], self.field)
                positions = locator_poly.roots()

                print('Positions: ', positions)

                indexes = []

                for p in positions:
                    i = -1
                    temp = p
                    while temp != 1:
                        temp = temp/self.field.primitive_element
                        i += 1
                    
                    indexes.append(i)

                print([int(s) for s in syndromes])
                syndrome_poly = galois.Poly(syndromes[::-1], self.field)

                deriv_coeffs = locator_poly.coeffs.copy()

                for i, coeff in enumerate(deriv_coeffs):
                    if (len(deriv_coeffs)-i-1) % 2 == 0:
                        deriv_coeffs[i] = 0


                locator_derivative = galois.Poly([0] + [int(coeff)  for coeff in deriv_coeffs[:-1]], self.field)
                print('Locator ', locator_poly)
                print('Derivative ', locator_derivative)

                print(syndrome_poly)

                evaluator = (locator_poly*syndrome_poly) % galois.Poly([1] + [0]*(2*self.t), self.field)

                print(evaluator)                    
                    #data[i] = 1 if data[i] == 0 else 0

                #print(self.field.)
                estimated_error_coeffs = [0] * self.n

                for count, p in enumerate(positions):
                    value = (evaluator(np.reciprocal(p)) / locator_derivative(np.reciprocal(p)))
                    print('bin value: ', bin(int(value)))
                    estimated_error_coeffs[indexes[len(indexes)-count-1]] = int(value)
                #estimated_error_coeffs[0] = 

                print('Estimated error coeffs: ', estimated_error_coeffs)
                error_estimator = galois.Poly(estimated_error_coeffs, self.field)
                print('Error estimator poly: ', error_estimator)
                corrected = received_poly - error_estimator
                print('Corrected? ', corrected)

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

    rs7_3 = RS(7, 3)

    #msg = np.random.randint(2, size=9)

    msg =  [0,1,0,1,1,0,1,1,1]
    encoded = rs7_3.encode(msg)
    print(encoded)
    encoded[3] = 1 if encoded[3] == 0 else 0
    encoded[0] = 1 if encoded[0] == 0 else 0
    print(encoded)
    rs7_3.correct(encoded)
