import random
import math

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def box_muller() -> float:
       u1 = np.random.uniform()
       u2 = np.random.uniform()

       return np.sqrt(-2*np.log(u1))*np.sin(2*np.pi*u2)


print(box_muller())


