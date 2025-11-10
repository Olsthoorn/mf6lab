# Steady-state wide-ditch head and seepage computation
#
#
# Analytic solution for flow in shallow aquifer on top of semi-confined
# regional aquifer with given head phi_deep.
# Situation
#
#                               |||||||||||| N2 |||||||||||||||||||
#     wide ditch ,h=hLR      |  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv |
#                            |                                      |
#    |~~~~~~~~  c_sb~~~~~~~~~~~                                      |
#    |          kD1                         kD2                     |
#    |//////////c1//////////////////////////c2///////////////////// |
#    |          phi_deep                         phi_deep                     |
#    |^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ |
#    |||||||||||||q1|||||||||||||||||||||||  q2 ||||||||||||||||||| |
#    |==== x1 ===>                                   <== x2 ========|
#    |<-------- b1 --------->|<-------------- b2 ------------------>|
#
# Ditch can be on either side, like the recharge.
#
# Flux through confining bed positive if upward (Dutch "kwel").
#
#
# TO 2010-12-27 2011-01-22 Matlab, 2025-07-08 Python
# %%
import numpy as np
import matplotlib.pyplot as plt

# %% Initialize solution specific

qsoll = -0.003         # [m/d] desired upward flux over entire width b1+b2
phi_deep = -1.0             # [ m ] prescribed head in regional aquifer
q = 0.001
w = 1.0


# First value left compartment, second value right compartment
N   = (0.001, 0.001)     # [m/d] reharge on the two compartments
hLR = (0.5,  0.5)       # [ m ] head in ditch
kD  = (15,  15)          # [m2/d] transmissivity of cover layer
C   = (20,  20)        # [ d ] resistance at bottom of cover layer
b   = (100, 100)        # [ m ] width of the two compartments
c_sb = (np.inf,  np.inf)     # [ m ] ditch bottom resistance, use Inf if no ditch
Q = 0.1 # [ m2/d ] discharge throug ditch fase.

x = (
    np.linspace(0, b[0], 100),  # x1 in first compartment
    np.linspace(0, b[1], 100)   # x2 in second
)

X = np.hstack((x[0] - x[0][-1], x[1][::-1]))

c = (
    1 / (1 / c_sb[0] + 1 / C[0]),
    1 / (1 / c_sb[1] + 1 / C[1])  # [ d ] effective resistance of compartment
)

lam =(
    np.sqrt(kD[0] * c[0]),
    np.sqrt(kD[1] * c[1])
)

# Representative heads in both compartments
phi = (
    (hLR[0] / c_sb[0] + phi_deep / C[0]) / (1 / c_sb[0] + 1 / C[0]),
    (hLR[1] / c_sb[1] + phi_deep / C[1]) / (1 / c_sb[1] + 1 / C[1])
)

def TshL(T, L):
    return T / L * np.tanh(L / T)

B = (
    np.cosh(b[0] / lam[0]) + TshL(kD[0], lam[0]) / TshL(kD[1], lam[1]) * np.cosh(b[1] / lam[1]),
    np.cosh(b[1] / lam[1]) + TshL(kD[1], lam[1]) / TshL(kD[0], lam[0]) * np.cosh(b[0] / lam[0])
)

beta = (
    lam[0] / b[0] * np.sinh(b[0] / lam[0]) / B[0],
    lam[1] / b[1] * np.sinh(b[1] / lam[1]) / B[1]
)

hb = (
    (
        (phi[0] + N[0] * c[0] + ((N[1] * c[1] + phi[1]) - (N[0] * c[0] + phi[0])) / B[0] * np.cosh(b[0] / lam[0]) +
         hLR[0] / (w * B[0]) * np.cosh(b[0] / lam[0]) / TshL(kD[1], lam[1])) /
        (1 + 1  / (w * B[0]) * np.cosh(b[0] / lam[0]) / TshL(kD[1], lam[1]))
    ),
    (
        (phi[1] + N[1] * c[1] + ((N[0] * c[0] + phi[0]) - (N[1] * c[1] + phi[1])) / B[1] * np.cosh(b[1] / lam[1]) + 
         hLR[1] / (w * B[1]) * np.cosh(b[1] / lam[1]) / TshL(kD[0], lam[0])) /
        (1 + 1  / (w * B[1]) * np.cosh(b[1] / lam[1]) / TshL(kD[0], lam[0])))
)

Qb = (
 (hb[0] - hLR[0]) / w,
 (hb[1] - hLR[1]) / w
)

s = (
    N[0] * c[0] + (N[1] * c[1] + phi[1] -
            (N[0] * c[0] + phi[0]) * Q / TshL(kD[1], lam[1]) / B[0] * np.cosh(x[0] / lam[0])) / B[0],
    N[1] * c[1] + (N[0] * c[0] + phi[0] -
            (N[1] * c[1] + phi[1]) * Q / TshL(kD[0], lam[0]) / B[1] * np.cosh(x[1] / lam[1])) / B[1]
)
   
savg = (
    N[0] * c[0] + beta[0] * (N[1] * c[1] + phi[1] -
            (N[0] * c[0] + phi[0]) - Q / TshL(kD[1], lam[1])),
    N[1] * c[1] + beta[1] * (N[0] * c[0] + phi[0] -
            (N[1] * c[1] + phi[1]) - Q / TshL(kD[0], lam[0]))
)

Savg = (b[0] * savg[0] + b[1] * savg[1]) / (b[0] + b[1])  # [ m2/d ] section-averaged seepage

h = (
    phi[0] + N[0] * c[0] + (savg[0] - N[0] * c[0]) * np.cosh(x[0] / lam[0]) / (lam[0] / b[0] * np.sinh(b[0] / lam[0])),
    phi[1] + N[1] * c[1] + (savg[1] - N[1] * c[1]) * np.cosh(x[1] / lam[1]) / (lam[1] / b[1] * np.sinh(b[1] / lam[1]))
)

dphi_dphi_deep = (
    c_sb[0] / (C[0] + c_sb[0]),
    c_sb[1] / (C[1] + c_sb[1])
)

dQ_dphi = (
    (dphi_dphi_deep[0] + (dphi_dphi_deep[1] - dphi_dphi_deep[0]) * np.cosh(b[0] / lam[0]) / B[0]) /
    (w + 1 / B[0]) * np.cosh(b[0] / lam[0]) / TshL(kD[1], lam[1]),
    (dphi_dphi_deep[1] + (dphi_dphi_deep[0] - dphi_dphi_deep[1]) * np.cosh(b[1] / lam[1]) / B[1]) / 
    (w + 1 / B[1]) * np.cosh(b[1] / lam[1]) / TshL(kD[0], lam[0])     
)

# phi_plus = (
    # (q - q[0]) / dQ_dphi[0] + phi_deep,
    # (q - q[0]) / dQ_dphi[1] + phi_deep
# )

# %% Coordinats to compute heads

x1=np.arange(0, b[0]) # x in first compartment from left to right
x2=np.arange(0, b[1]) # x in second compartment from right to left

x =np.unique([x1, b[1]- x2])  # true real world x along cross section from center left

fig, ax = plt.subplots(figsize=(10, 6))
ax.grid(True)
ax.set_xlabel('x [m]')
ax.set_ylabel('head [m]')
ax.set_title(f'wide ditch, q={q:.3f} m/d')

ax.plot(-x1, h[0],'r', label='h1') 
ax.plot(+x2, h[1],'r', label='h2')
ax.plot(-x1, savg[0] * np.ones_like(x1),'b', label='savg[0]')
ax.plot(+x2, savg[1] * np.ones_like(x2),'b', label='savg[1]')

ax.plot(-x1, phi_deep * np.ones_like(x1), 'k')
ax.plot(+x2, phi_deep * np.ones_like(x2), 'k', label='phi_deep')
ax.legend()


# %% Computation of heads and mean heads and loop once to match q with
#  desired section-averaged seepage rate



# NStep=2  # We need exactly two steps to get the desired flux
# hmean=np.zeros((1, 2)) * np.nan
# Out = np.zeros((NStep,2)) * np.nan # [ m ] output array for phi_deep and q
# Out[0,:] = (phi_deep, 0)  # initial values

# for i in range(1, NStep):

#     phi = (hLR/c_sb + phi_deep/C) / (1/c_sb + 1/C)  # net fixed head in compartments

#     # head as function of x1 and x2
#     h1=N[0] * c[0] + ((N[1] * c[1] + phi[1]) - (N[0] * c[0] + phi[0])) * np.cosh(x1 / lam[1]) / B[1] + phi[0]
#     h2=N[1] * c[1] + ((N[0] * c[0] + phi[0]) - (N[1] * c[1] + phi[1])) * np.cosh(x2 / lam[2]) / B[2] + phi[1]

#     # average head in both compartments
#     hmean[0]=N[0] * c[0] + ((N[1] * c[1] + phi[1]) - (N[0] * c[0] + phi[0])) * (lam[1] / b[0]) * np.sinh(b[0] / lam[1]) / B[1] + phi[0]
#     hmean[1]=N[1] * c[1] + ((N[0] * c[0] + phi[0]) - (N[1] * c[1] + phi[1])) * (lam[2] / b[1]) * np.sinh(b[1] / lam[2]) / B[2] + phi[1]

#     # section-averaged seepage
#     q=b[0] / np.sum(b) * (phi_deep - hmean[0]) / C[1] + b[1] / np.sum(b) * (phi_deep - hmean[1]) / C[2] # seepage

#     # show
#     print(f'phi_deep({i})={phi_deep:10f}  q={q:10g}')

#     # keep
#     Out[i,:] = [phi_deep, q]

#     # update phi_deep for next loop
#     phi_deep += (qsoll - q) / dqdphi

# phi_deep=Out[-1, 0] # set phi_deep back to last used value

# phi_deep=Out[-1, 0] # set phi_deep back to last used value

# %% Plot results

