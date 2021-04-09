import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
import matplotlib.animation as ani

# --------------------------------- #
NBATH = 3
NSITES = NBATH + 1 # 1 QD site
U = 4 # local coulomb interaction
e = 1 # electron charge
E = -2 # QD energy level
Ek = np.linspace(-NBATH, NBATH, NBATH) # evenly spaced Bath levels around fermi level = 0
Vk = [1]*(NBATH+1) # Bath-Lead hoppings
# --------------------------------- #