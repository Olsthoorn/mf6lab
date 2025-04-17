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

# %% Grid setup
nx, nz = 1000, 100
x = np.linspace(0, 1000, nx)
z = np.linspace(0, 100, nz)
grid_x, grid_z = np.meshgrid(x, z)
gr = Grid(x, None, z, axial=False)
extent = (gr.x[0], gr.x[-1], gr.z[-1], gr.z[0])

# %% Variogram parameters
name = 'gaus'
name = 'wim'
KG = 0.89 # m/d
mu = np.log(KG)
var_ln_k = 6.6          # variance (sill)
# var_ln_k = 1.5          # variance (sill) homogeneous

len_x = 10.2       # correlation length in x (range)
len_z = 1.5        # correlation length in y (range)
model = gs.Exponential(dim=2, var=var_ln_k, len_scale=[len_x, len_z])

k_field_pars = {'name': name,
                'kG': KG,
                'var_ln_k':var_ln_k,
                'Ih':len_x,  # Horizontal integration scale [m]
                'Iv':len_z, # Vertical integration scale [m]
                'theta': 0.31,
                'J': 0.0036 # mean gradient
              }

kfp = k_field_pars

k_field_str = fr"$k_G={kfp['kG']:.3g}\,m/d,\,var(\ln(k)={kfp['var_ln_k']:.3g},\,I_h={kfp['Ih']:.3g}\,m,\,I_v={kfp['Iv']:.3g}\,m,\,\theta={kfp['theta']:.2g},\,J={kfp['J']:.2g}$"

kfp['k_field_str'] = k_field_str

# %%  Generate field
srf = gs.SRF(model, seed=42)
field_ln_k = srf((gr.xm, gr.zm), mesh_type='structured').T
mim_ln_k = np.zeros_like(field_ln_k)

# %% Fill the block_k_field

if name == 'wim':
    k1, k2, mu = 5., 500., 0.
    nx_block, nz_block, nx_incl, nz_incl= 100, 5, 20, 1
    mim_ln_k[:, :] = np.log(k1)    
    for ixB in range(0, gr.nx, nx_block):
        for izB in range(0, gr.nz, nz_block):
            ixIncl = np.random.choice(nx_block - nx_incl)
            izIncl = np.random.choice(nz_block - nz_incl)
            i1, i2 = ixB + ixIncl, ixB + ixIncl + nx_incl
            j1, j2 = izB + izIncl, izB + izIncl + nz_incl          
            mim_ln_k[j1:j2, i1:i2] = np.log(k2)
    k_field_str = fr"$k1={k1}\,m/d,\,k_2={k2}\,m/d$, nx,nz blocks=({nx_block},{nz_block}), nx, nz inclusions=({nx_incl},{nz_incl})"
    kfp['k_field_str'] = k_field_str    
elif name != 'gaus':
    for iz in range(gr.nz):
        X = -np.random.rand(1) * 2 * len_x + np.arange(0, x[-1] + 2 * len_x, 2 * len_x)
        Idx = np.searchsorted(gr.xm, X)

        for i1, i2 in zip(Idx[:-1], Idx[1:]):
            im = (i1 + i2) //2
            if im >= nx - 1:
                im = -1
            if i2 >= nx - 1:
                i2 = None   
            mim_ln_k[iz, i1:i2] = field_ln_k[iz, im]
else:
    mim_ln_k = field_ln_k

mim_ln_k += mu

k_field = np.exp(mim_ln_k)

# %% Plot
if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(10, 4))    
    mappable = ax.imshow(mim_ln_k, extent=extent, aspect=1, origin='lower', cmap='viridis')
    c = plt.colorbar(mappable,label="Log-conductivity", orientation='horizontal', ax=ax)
    c.set_label("ln(k)")
    ax.set(title="2D Random ln(k) field\n" + k_field_str,
           xlabel='x [m]', ylabel='y [m]')
    plt.tight_layout()
    plt.show()

# %%
