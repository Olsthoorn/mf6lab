lrcMilk  = gr.lrc_recarray(*np.array(pr['milkInjPnt']).T) # global coords milk injection point

IDOMAIN.ravel()[gr.Iglob(lrcMilk)] = pr['iMlkInjPnt'] # mark milk injection cells
