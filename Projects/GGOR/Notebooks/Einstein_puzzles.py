# %% [markdown] Einstein puzzles
#
# Einstein puzzels zijn logische puzzels waarin je op basis
# van in zinnen verpakte gegevens unieke combinaties
# moet bepalen van een aantal categoriëen.
#
# Vragen:
# 
# 1. Anna heeft niet de hond.
# 2. Wie de kat heeft, graag een groene trui.
# 3. Bram draagt geen blauwe trui.
# 4. Clara heeft de vis.
# %%

import numpy as np

# %%

dtype = [('pers', str ), ('dier', str), ('kleur', str)]
opl = np.zeros((3, 3), dtype=dtype)

# %%
personen = ['Anna', 'Bram', 'Clara']
dieren = ['kat', 'hond', 'vis']
kleuren = ['blauw', 'groen', 'rood']

# %%
for wie in personen:
    for dier in dieren:
        for kleur in kleuren:
            for opl in oplossingen:
                if opl['wie'] is None and opl['dier'] is None:
                    opl['wie'] = 'Anna'
                elif not(wie == 'Anna' and not(dier == 'hond')):
                
            #print(k)
            if opl['w'] is None and opl['d'] != 'hond':
                opl['w'] = 'Anna'
                continue
            if opl['d'] == 'kat':
                opl['k'] = 'groen'
                continue
            if opl['w']== 'Bram':
                and (not (kleur == 'blauw'))):
                break
            if not ((wie == 'Clara') and (dier == 'vis')):
                break
            print(wie, dier, kleur)
# %%
