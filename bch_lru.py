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


# if __name__ == '__main__':
#     bch15_7 = BCH(4, galois.Poly([1, 0, 0, 1, 1], galois.GF(2)), [1,1,1,0,1,0,0,0,1], 2)

#     true = 0
#     total = 0

#     for i in range(0,15):
#         for j in range(0,15):

#             if i == j:
#                 continue

#             total += 1


#             msg = np.random.randint(2, size=7)
#             codeword = bch15_7.encode(msg)
#             modified = codeword.copy()
#             print(codeword)
#             modified[i] = 1 if modified[i] == 0 else 0
#             modified[j] = 1 if modified[j] == 0 else 0
#             #codeword[30] = 1 if codeword[30] == 0 else 0

#             corrected = bch15_7.correct(modified)
#             #print(modified)
#             #print(len(bch127.correct(modified)))
#             decoded = bch15_7.decode(corrected)
#             #print(len(decoded))
#             #decoded[0] = 1 if decoded[0] == 0 else 0
#             #print(decoded == msg)
#             if all(decoded == msg):
#                 true += 1

#     print(true, total)

from utils import chunks
from box_muller import box_muller
import sys
import time
import csv
import json
import os

if __name__ == '__main__':

    simulation_data = {}
        
    bch_code = bch15_7()

    BCH_BITS = int(1000000 - 1000000%7)

    snr_in_db_range = np.arange(0, 9, 0.5)

    bch_bers = []

    start = time.time()

    for index in range(len(snr_in_db_range)):

            local_bers = []

            snr_in_db = snr_in_db_range[index]

            simulation_data[snr_in_db] = {}

            snr_linear = 10**(snr_in_db/10)

            noise_spectral_dens = 1/snr_linear

            for i in range(10):
                simulation_data[snr_in_db][i] = {}

                uncoded = np.random.randint(0, 2, BCH_BITS, int)

                simulation_data[snr_in_db][i]['uncoded_message'] = uncoded.tolist()

                msg = []

                for chunk in chunks(uncoded, 7):
                        msg += bch_code.encode(tuple(chunk))

                simulation_data[snr_in_db][i]['encoded_message'] = msg

                signal = 1 - 2*np.asarray(msg)

                print('n0', noise_spectral_dens)

                randoms = [box_muller() for _ in range(len(msg))]

                noise = np.sqrt(noise_spectral_dens/2)*np.asarray(randoms)

                received = signal + noise

                demodulated = [1 if r < 0 else 0 for r in received]

                simulation_data[snr_in_db][i]['demodulated_message'] = demodulated
                
                decoded = []

                for dem_chunk in chunks(demodulated, 15):
                        decoded += bch_code.decode(tuple(bch_code.correct(tuple(dem_chunk))))

                simulation_data[snr_in_db][i]['decoded_message'] = decoded

                n_errors = sum(decoded != uncoded)
                ber = n_errors/BCH_BITS

                simulation_data[snr_in_db][i]['number_of_errors'] = int(n_errors)
                simulation_data[snr_in_db][i]['bit_error_rate'] = float(ber)

                local_bers.append(ber)


            bch_bers.append(np.average(local_bers))
    end = time.time()
    print('Elapsed: ', end - start)

    if not os.path.exists('bch'):
        os.makedirs('bch')

    with open(f'bch/bch_aggregated_bers_{end}.csv', mode='w') as file:
            writer = csv.writer(file)
            writer.writerows(map(lambda x: [x], bch_bers))

    with open(f'bch/bch_raw_data_{end}.json', mode='w') as file:
        json.dump(simulation_data, file)
