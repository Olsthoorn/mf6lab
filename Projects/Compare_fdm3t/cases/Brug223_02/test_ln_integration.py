# See if integrating along a log axis instead a  linear axis is done
# correctly. The conclusion is, that that is indeed the case. Also, that
# the quad function does not work well if the axis runs over a buch of
# log cycles. In that case integration should  be formulate along the
# log axis directly, turning it into a linear axis.
#
# exp(-s x) is used as a sample function
#

# TO 2024-12-28

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def fu(u, s):
    return np.exp(-s * u)

def fz(z, s):
    u = np.exp(z)
    return u * np.exp(-s * u)

def simpson(f, x, args):
    F = f(x, *args)
    return np.sum(np.diff(x) * 0.5 * (F[:-1] + F[1:]))

u = np.logspace(-3, 6,1000)
z = np.log(u)

u0, u1 = u[[0, -1]]
z0, z1 = z[[0, -1]]

s = 4.0

fu1, efu1 = quad(fu, u0, u1, (s))
fz1, efz1 = quad(fz, z0, z1, (s))

print(f"fu1 = {fu1}, efu1 = {efu1}")
print(f"fz1,= {fz1}, efz1 = {efz1}")
print(f"fu(simpson) = {simpson(fu, u, (s,))}")
print(f"fz(simpson) = {simpson(fz, z, (s,))}")
print('hello world')

# The answer proves
# 1 that the procedure is correct
# 2 quad can't handle the linear scale when dealing with a number of log cycles.