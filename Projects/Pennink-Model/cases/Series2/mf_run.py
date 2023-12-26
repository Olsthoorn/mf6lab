# Run Modflow etc for the current Case
import os
import logging
import matplotlib.pyplot as plt
import numpy as np

cws = '/Users/Theo/GRWMODELS/python/Pennink-Model/cases/Series2'
os.chdir(cws)

import mf_adapt
from mf6lab import mf_setup

logging.info("Running from {}".format(cws))

gr = mf_adapt.gr

fp_packages, model_dict, use_models, use_packages = mf_setup.mf_setup()

sim = fp_packages.get('Simsim')

# Write simulation
sim.write_simulation(silent=False)

# Run simulation
success, buff = sim.run_simulation(silent=False)

print('Running success = {}'.format(success))
if not success:
    print(buff)
    logging.critical("Buffer printed because MODFLOW did not terminate normally.")
    raise Exception('MODFLOW did not terminate normally.')

# TODO Indicate the seepage zone.

print("Done")