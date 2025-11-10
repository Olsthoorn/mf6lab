# Fiori at al kfield Geneateor

# All bocks have the same size (2R)
# The blocks in adjacent rows are somewhat shifted relative to allow paprticles
# to naturally step into adjacent rows.
# There is a autocorrelation structure defined by the integration scales.
# We use the lognormal conductiviy distribution together with the
# integration scale to sample condictivities.

# We can fill row by row, from below using the horizontal and vertical integration
# scales, i.e. by sampling from the k-distribution with a reduced sigma because
# the distance to the adjacent cell is smaller than the intgration scale.

# %% imports
import numpy as np
import matplotlib.pyplot as plt
import gstools as gs
from fdm.mfgrid import Grid
import random

# %% Grid setup
nxp, nzp = 2002, 101
x = np.linspace(-0.5, 2000.5, nxp)
z = np.linspace(0, 50, nzp)
grid_x, grid_z = np.meshgrid(x, z)
gr = Grid(x, None, z, axial=False)
extent = (gr.xm[0], gr.xm[-1], gr.z[-1], gr.z[0])

kG, var_ln_k, Ih, Iv, theta, J = 0.89, 6.6, 10.2, 1.5, 0.31, 0.0036

# %% Variogram parameters
case_name, fJ = 'exp1',   1.000
# case_name, fJ = 'exp2',   0.708
# case_name, fJ = 'exp4',   0.494
# case_name, fJ = 'exp8',   0.308
# case_name, fJ = 'exp16',  0.224
# case_name, fJ = 'exp32',  0.175
# case_name, fJ = 'exp64',  0.126
# case_name, fJ = 'exp128', 0.115

if case_name.startswith('exp'):
    Iv = 1.0
    Ih = float(case_name[3:])

k_field_pars = {'name': case_name,
                'kG': kG,
                'var_ln_k':var_ln_k,
                'Ih': Ih,  # Horizontal integration scale [m]
                'Iv': Iv, # Vertical integration scale [m]
                'theta': theta,
                'J': J * fJ # mean gradient
              }

kfp = k_field_pars

# Basis k_field_str:
k_field_str = fr"$k_G={kfp['kG']:.3g}\,m/d,\,var(\ln(k)={kfp['var_ln_k']:.3g},\,I_h={kfp['Ih']:.3g}\,m,\,I_v={kfp['Iv']:.3g}\,m,\,\theta={kfp['theta']:.2g},\,J={kfp['J']:.2g}$"

kfp['k_field_str'] = k_field_str

# %% Fill the block_k_field
if case_name.startswith('exp'):
    # Exponential random field with integration lengths Ih and Iv and var_ln(k) given.
    model = gs.Exponential(dim=2, var=var_ln_k, len_scale=[Ih, Iv])
    srf = gs.SRF(model, seed=42)
    field_ln_k = srf((gr.xm, gr.zm), mesh_type='structured').T
    field_ln_k += np.log(kG)
    k_field = np.exp(field_ln_k)
    kfp['k_field_str'] = k_field_str

elif case_name == 'random': # k is randomly selected from lognormal distribution
    # This case aims to show the breakthrough curves and statistics for a
    # purely random conductivity field, where the k in each cells is randomly
    # selected from a lognormal distribution with given kG and sigma.
    mu = np.log(kG)
    sigma = np.sqrt(var_ln_k)
    # k_field = np.exp(np.random.normal(mu, sigma, size=gr.shape))[:, 0, :]
    k_field = np.random.lognormal(mu, sigma, size=gr.shape)[:, 0, :]
    
elif case_name.startswith('mim'): # Multi Indicator Model (Fiori at al 2013)
    # This aims to follow the setup by Fiori et al (2013) in thir Multi Indicator Model.
    # The cross section consists of blocks of length 2Ih with conductivitiy values
    # drawn from the lognorml exponential random field with given properties.
    # The blocks are shifted with respect to each other.
    model = gs.Exponential(dim=2, var=var_ln_k, len_scale=[Ih, Iv])
    srf = gs.SRF(model, seed=42)
    field_ln_k = srf((gr.xm, gr.zm), mesh_type='structured').T
    mim_ln_k = np.zeros_like(field_ln_k)

    for iz in range(gr.nz):
        # Horizontally shift the blocks in a layer by choosing a random offset for this layer
        X = -np.random.rand(1) * 2 * Ih + np.arange(0, x[-1] + 2 * Ih, 2 * Ih)
        Idx = np.searchsorted(gr.xm, X)

        # Then for each block, select the  k value from the Gauss ln(k)-field.
        # The value chosen is the gauss field value at the center of the block.
        for i1, i2 in zip(Idx[:-1], Idx[1:]):
            im = (i1 + i2) //2
            if im >= gr.nx:
                im = -1
            if i2 >= gr.nx:
                i2 = None   
            mim_ln_k[iz, i1:i2] = field_ln_k[iz, im]
    mim_ln_k += np.log(kG)
    k_field = np.exp(mim_ln_k)
    kfp['k_field_str'] = k_field_str

elif case_name == 'WdL': # Wim de Lange's conceptual domains model
    # Wim's aquifer with background kG and k1 and k2 as high and low k inclusions.
    # Ih = 10.2,  Iv=1.5, La = 6 * Ih = nx_block, and Da = 6 * Iv = nz_block.
    # such that the inclusion fits entirely within the domain.
    sigma = np.sqrt(var_ln_k)
    k1, k2 = kG * np.exp(sigma), kG /np.exp(sigma)
    k_field = gr.const(kG)[:, 0, :]
    nx_block, nz_block = int(6 * Ih), int(6 * Iv)
    nx_incl, nz_incl= int(2 * Ih), int(2 * Iv)
    for ixB in range(0, gr.nx + 1, nx_block):
        for izB in range(0, gr.nz + 1, nz_block):
            ix1, ix2 = random.sample(range(nx_block - nx_incl), 2)
            iz1, iz2 = random.sample(range(nz_block - nz_incl), 2)
            
            i1, i2 = ixB + ix1, ixB + ix1 + nx_incl
            j1, j2 = izB + iz1, izB + iz1 + nz_incl          
            k_field[j1:j2, i1:i2] = k1
            
            i1, i2 = ixB + ix2, ixB + ix2 + nx_incl
            j1, j2 = izB + iz2, izB + iz2 + nz_incl
            k_field[j1:j2, i1:i2] = k2         

    k_field_str = fr"$k_G={kG},\,m/d,\,k1={k1:.3g}\,m/d,\,k_2={k2:.3g}\,m/d$, nx,nz blocks=({nx_block},{nz_block}), nx, nz inclusions=({nx_incl},{nz_incl})"
    kfp['k_field_str'] = k_field_str

else:
    raise ValueError(f"Unknown case {case_name}")

kfp['k_field'] = k_field

# %% Plot
if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(10, 4))    
    mappable = ax.imshow(np.log(k_field), extent=extent, aspect=1, origin='lower', cmap='viridis')
    c = plt.colorbar(mappable,label="Log-conductivity", orientation='horizontal', ax=ax)
    c.set_label("ln(k)")
    ax.set(title="2D Random ln(k) field\n" + k_field_str,
           xlabel='x [m]', ylabel='y [m]')
    plt.tight_layout()
    plt.show()

# %%
