import numpy as np
from scipy.special import gamma

s = 1
for i in range(365, 342, -1):
    s *= i / 365
    print(i, s)
    
print()
potenz = (353 / 365) ** 22
print('353^22 = ', potenz)

print(potenz / s)
