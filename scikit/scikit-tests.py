from sk_dsp_comm import fec_block, fec_conv
from numpy import random


msg = random.randint(0, 2, 4)
print(msg)

hamming7_4 = fec_block.fec_hamming(3)

encoded = hamming7_4.hamm_encoder(msg)
print(encoded)
encoded[1] = 1 if encoded[1] == 0 else 0
print(encoded)
decoded = hamming7_4.hamm_decoder(encoded.astype(int))
print(decoded)