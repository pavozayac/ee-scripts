from math import log2
from galois.field.conway import conway_poly
from galois.field.poly_functions import primitive_element
import numpy as np 
import galois
from numpy.core.defchararray import encode, index
from numpy.core.function_base import logspace
from numpy.lib.polynomial import polydiv

class RS():
    
    def __init__(self, n: int, k: int) -> None:
        self.n = n
        self.k = k
        self.order = int(log2(n+1))
        self.t = int((n-k)/2)


        self.field = galois.GF(2**self.order, conway_poly(2, self.order))
        self.field.display('power')

        a = self.field.primitive_element
        self.generator = galois.Poly.Roots(roots=[a**i for i in range(1, 2*self.t+1)], field=self.field)
        print(self.generator.roots())
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

        shifted = msg_poly * galois.Poly([1, *[0]*(self.n-self.k)], self.field)

        parity_poly = shifted % self.generator

        codeword = parity_poly + shifted
        print('Codeword: ', codeword)
        print('Codeword as binary: ', [bin(coeff) for coeff in codeword.coeffs])

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
        syndrome_poly = received_poly % self.generator

        syndromes = []

        for i in range(1, 2*self.t+1):
            syndromes.append(syndrome_poly(self.field.primitive_element**i))

        print('syndromes: ', syndromes)
        for v in range(self.t, 0, -1):
            print("\n V: ", v)
            m = []

            #galois.log.log_naive

            for i in range(v):
                m.append(syndromes[i:v+i])
                
            if np.linalg.det(self.field(m)) != 0:
                svs = self.field([[x] for x in syndromes[v:2*v]])

                print(self.field(svs))
                print(m)
                print(np.linalg.inv(self.field(m)))

                loc_coeffs = np.linalg.inv(self.field(m)) @ self.field(svs)
                
                print('loc coeffs ', loc_coeffs)

                locator_poly = galois.Poly([*loc_coeffs.flatten(), 1], self.field, order='asc')

                print('loc poly: ', locator_poly)
                positions = locator_poly.roots()    

                print('Positions: ', positions)

                indexes = []

                for p in positions:
                    indexes.append(np.log(p**-1))

                print('Indexes ', indexes)

                pos_matrix = []

                for i in range(v):
                    pos_matrix.append([pos**(i+1) for pos in positions])

                print('Position matrix ', pos_matrix)

                values = np.linalg.inv(self.field(pos_matrix)) @ self.field([[syn] for syn in syndromes[0:v]])

                print('Values: ', values)
                print('Values flattenned: ', values.flatten())

                error_coeffs = [0] * self.n

                for count, i in enumerate(indexes):
                    error_coeffs[i-1] = values.flatten()[count]

                error_poly = galois.Poly(error_coeffs, self.field)

                print('Error poly: ', error_poly)

                corrected = received_poly+error_poly

                print('Corrected: ', corrected)
                """
                
                alpha = self.field.primitive_element

                print([int(s) for s in syndromes])

                deriv_coeffs = locator_poly.coeffs.copy()

                for i in range(len(deriv_coeffs)):
                    deriv_coeffs[i] = deriv_coeffs[i]*i


                locator_derivative = galois.Poly([0] + [int(coeff)  for coeff in deriv_coeffs[:-1]], self.field)
                print('Locator ', locator_poly)
                print('Derivative ', locator_derivative)

                print('Syndrome poly', syndrome_poly)

                #evaluator = (locator_poly*syndrome_poly) % galois.Poly.Degrees([2*self.t], field=self.field)

                eval_coeffs = [0] * (2*self.t+1)

                for i in range(2*self.t+1):
                    eval_coeffs[i] = locator_poly(alpha**i)*syndrome_poly(alpha**i)

                print('evaluator coeffs: ', eval_coeffs)

                evaluator = (locator_poly*syndrome_poly) % galois.Poly.Degrees([2*self.t], field=self.field)
                eval_alt = galois.Poly(eval_coeffs, field=self.field, order='asc')

                print('evaluator', evaluator)
                print('evaluator alt', eval_alt)
                #for i in range(1, v+1):
                #    syn_vec = [1] + syndromes[:i-1]
                #    loc_vec = [[]] + [[1]]

                #print('Evaluator', evaluator)                    
                    #data[i] = 1 if data[i] == 0 else 0

               # X_mat = []

                #for i in range(1, v+1):
                #    X_mat.append([pos**i for pos in positions])

                
                #for count, row in enumerate(X_mat):
                #   X_mat[count].append(syndromes[count])

                # print("x_mat", X_mat)

                #dec_a, u = galois.array.lu_decompose(X_mat)
                #mags = np.linalg.inv(self.field(X_mat)) @ self.field([[syn] for syn in syndromes[:v][::-1]])

                #print('mags from solving', mags)
                #print(self.field.)
                estimated_error_coeffs = [0] * self.n

                
                print('recipr', positions[0], positions[0]**-1, positions[0]*positions[0]**-1)

                for i in indexes:
                    print('inner loc_der eval',  locator_derivative(alpha**-i))
                    value = evaluator(alpha**-i) / locator_derivative(alpha**-i)
                    value_alt =  (eval_alt(alpha**-i) / locator_derivative(alpha**-i))
                    print('bin value: ', bin(int(value)))
                    print('bin value alt: ', bin(int(value_alt)), value_alt, value_alt * (alpha**i))
                    estimated_error_coeffs[i-1] = int(value)
                #estimated_error_coeffs[0] = 

                print('Estimated error coeffs: ', estimated_error_coeffs)
                error_estimator = galois.Poly(estimated_error_coeffs, self.field, order='desc')
                print(error_estimator.roots())
                print('Error estimator poly: ', error_estimator)
                corrected = received_poly + error_estimator
                print('Corrected? ', corrected)
                """
            

                

                corrected_coeff_binstr = [bin(int(coeff))[2:] for coeff in corrected.coeffs]

                for index in range(len(corrected_coeff_binstr)):
                    corrected_coeff_binstr[index] = [int(bit) for bit in corrected_coeff_binstr[index]]
                    corrected_coeff_binstr[index] = (self.order-len(corrected_coeff_binstr[index]))*[0] + corrected_coeff_binstr[index]
                
                return [bit for sym in corrected_coeff_binstr for bit in sym]

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
    original_encoded = rs7_3.encode(msg)
    encoded = original_encoded.copy()
    print('Encoded\n', encoded)
    #encoded[5] = 1 if encoded[5] == 0 else 0
    encoded[2] = 1 if encoded[2] == 0 else 0
    #encoded[8] = 1 if encoded[8] == 0 else 0
    #encoded[18] = 1 if encoded[18] == 0 else 0

    print('Encoded with errors\n', encoded)

    corrected = rs7_3.correct(encoded)

    print(corrected)
    print(corrected == original_encoded)
