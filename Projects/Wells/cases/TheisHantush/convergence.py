# Read and analyze convergence
import os
import logging
import pandas as pd

from settings import dirs

os.chdir(dirs.case) #os.chdir(os.path.dirname(__file__))
logging.info("Running from {}".format(os.getcwd()))

inner = pd.read_csv(os.path.join(dirs.SIM, 'inner.csv'), header=0)
outer = pd.read_csv(os.path.join(dirs.SIM, 'outer.csv'), header=0)


