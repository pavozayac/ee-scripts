from math import log2
from galois import conway_poly
#from galois.field.poly_functions import primitive_element
import numpy as np 
import galois
from numpy.core.defchararray import decode, encode, index
from numpy.core.function_base import logspace
from numpy.lib.polynomial import polydiv
from functools import lru_cache

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
        #print(self.generator.roots())
        #print(self.generator)
        #t is the error-correcting capacity of the code

    @lru_cache(maxsize=None)
    def encode(self, data: tuple) -> list:
        data = list(data)

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
        #print('Codeword: ', codeword)
        #print('Codeword as binary: ', [bin(coeff) for coeff in codeword.coeffs])

        binstr_coeffs = [bin(int(coeff))[2:] for coeff in codeword.coeffs]

        for index in range(len(binstr_coeffs)):
            binstr_coeffs[index] = [int(bit) for bit in binstr_coeffs[index]]
            binstr_coeffs[index] = (self.order-len(binstr_coeffs[index]))*[0] + binstr_coeffs[index]

        enc = [bit for sym in binstr_coeffs for bit in sym]
        return [0]*(21-len(enc)) + enc 

        '''
        msg = galois.Poly(data, self.field)
        codeword: galois.Poly = self.generator * msg

        encoded = [int(b) for b in codeword.coeffs]
        padded = [0]*(127-len(encoded)) + encoded
        return padded

        '''
    @lru_cache(maxsize=None)
    def correct(self, data: tuple) -> list:
        data = list(data)

        symbols = []

        for i in range(int(len(data)/self.order)):
            symbols.append(data[i*self.order:i*self.order+self.order])

        coeffs = []

        for sym in symbols:
            coeffs.append(self.field.Elements()[int(''.join(str(bit) for bit in sym), 2)])

        #print(coeffs)
        received_poly = galois.Poly(coeffs, self.field)
        #print('Received: ', received_poly)

        #received = galois.Poly(data, self.field)
        syndrome_poly = received_poly % self.generator

        syndromes = []

        for i in range(1, 2*self.t+1):
            syndromes.append(syndrome_poly(self.field.primitive_element**i))

        #print('syndromes: ', syndromes)

        try: 
            for v in range(self.t, 0, -1):
                #print("\n V: ", v)
                m = []

                    
                for i in range(v):
                    m.append(syndromes[i:v+i])
                    
                if np.linalg.det(self.field(m)) != 0:
                    svs = self.field([[x] for x in syndromes[v:2*v]])

                    #print(self.field(svs))
                    #print(m)
                    #print(np.linalg.inv(self.field(m)))

                    loc_coeffs = np.linalg.inv(self.field(m)) @ self.field(svs)
                    
                    #print('loc coeffs ', loc_coeffs)

                    locator_poly = galois.Poly([*loc_coeffs.flatten(), 1], self.field, order='asc')

                    #print('loc poly: ', locator_poly)
                    positions = locator_poly.roots()    

                    #print('Positions: ', positions)

                    indexes = []

                    for p in positions:
                        indexes.append(np.log(p**-1))

                    #print('Indexes ', indexes)

                    pos_matrix = []

                    for i in range(v):
                        pos_matrix.append([pos**(i+1) for pos in positions])

                    #print('Position matrix ', pos_matrix)

                    values = np.linalg.inv(self.field(pos_matrix)) @ self.field([[syn] for syn in syndromes[0:v]])

                    #print('Values: ', values)
                    #print('Values flattenned: ', values.flatten())

                    error_coeffs = [0] * self.n

                    for count, i in enumerate(indexes):
                        error_coeffs[i-1] = values.flatten()[count]

                    error_poly = galois.Poly(error_coeffs, self.field)

                    #print('Error poly: ', error_poly)

                    corrected = received_poly+error_poly

                    #print('Corrected: ', corrected)

                    corrected_coeff_binstr = [bin(int(coeff))[2:] for coeff in corrected.coeffs]

                    for index in range(len(corrected_coeff_binstr)):
                        corrected_coeff_binstr[index] = [int(bit) for bit in corrected_coeff_binstr[index]]
                        corrected_coeff_binstr[index] = (self.order-len(corrected_coeff_binstr[index]))*[0] + corrected_coeff_binstr[index]
                    
                    enc = [bit for sym in corrected_coeff_binstr for bit in sym]
                    return [0]*(21-len(enc)) + enc

            return data
        except:
            return data

    def decode(self, data: list) -> list:        
        return data[0:9][::-1]
    

from utils import chunks
from box_muller import box_muller
import sys
import time
import csv
import json
import os

if __name__ == '__main__':

    simulation_data = {}

    rs_code = RS(7,3)

    RS_BITS = int(1000000 - 1000000%7)

    snr_in_db_range = np.arange(0, 9, 0.5)

    rs_bers = []

    start = time.time()

    for index in range(len(snr_in_db_range)):
        local_bers = []

        snr_in_db = snr_in_db_range[index]

        simulation_data[snr_in_db] = {}

        snr_linear = 10**(snr_in_db/10)

        noise_spectral_dens = 1/snr_linear

        for i in range(10):
            simulation_data[snr_in_db][i] = {}

            uncoded = np.random.randint(0, 2, RS_BITS, int)

            simulation_data[snr_in_db][i]['uncoded_message'] = uncoded.tolist()

            msg = []

            for chunk in chunks(uncoded, 9):
                    msg += rs_code.encode(tuple(chunk))

            simulation_data[snr_in_db][i]['encoded_message'] = msg

            signal = 1 - 2*np.asarray(msg)

            print('n0', noise_spectral_dens)

            randoms = [box_muller() for _ in range(len(msg))]

            noise = np.sqrt(noise_spectral_dens/2)*np.asarray(randoms)

            received = signal + noise

            demodulated = [1 if r < 0 else 0 for r in received]

            simulation_data[snr_in_db][i]['demodulated_message'] = demodulated
            
            decoded = []

            for dem_chunk in chunks(demodulated, 21):
                    decoded += rs_code.decode(tuple(rs_code.correct(tuple(dem_chunk))))

            simulation_data[snr_in_db][i]['decoded_message'] = decoded

            n_errors = sum(decoded != uncoded)
            ber = n_errors/RS_BITS

            simulation_data[snr_in_db][i]['number_of_errors'] = int(n_errors)
            simulation_data[snr_in_db][i]['bit_error_rate'] = float(ber)

            local_bers.append(ber)

        rs_bers.append(np.average(local_bers))

    print(rs_bers)
    end = time.time()
    print('Elapsed: ', end - start)

    if not os.path.exists('rs'):
        os.makedirs('rs')

    with open(f'rs/rs_aggregated_bers_{end}.csv', mode='w') as file:
            writer = csv.writer(file)
            writer.writerows(map(lambda x: [x], rs_bers))

    with open(f'rs/rs_raw_data_{end}.json', mode='w') as file:
        json.dump(simulation_data, file)


