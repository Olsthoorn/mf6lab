# Run Modflow etc for the current Case
import os
import logging

from src.mf_setup import mf_setup

logging.info("Running from {}".format(os.getcwd()))

fp_packages, model_dict, use_models, use_packages = mf_setup()

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
else:
    print("MF6 terminated normally, run mf_analyze to visualize the results.")