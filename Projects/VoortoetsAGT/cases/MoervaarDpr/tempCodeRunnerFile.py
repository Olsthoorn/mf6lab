# %% Get psi and suitable levels
P =  mf_adapt.rch * gr.dx.sum() / 4 # estimate of psi extremes
psi = gr.psi_row(fflows['frf'][iper], row=0)
levels = get_contour_levels(-P, P, 100)
