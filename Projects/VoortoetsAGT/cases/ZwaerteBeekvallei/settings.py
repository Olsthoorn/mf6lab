import os

from src.mf6tools import  Dirs

HOME = '/Users/Theo/GRWMODELS/python/mf6lab/Projects/VoortoetsAGT/'
assert os.path.isdir(HOME), "Can't find the directory {}".format(HOME)

LENGTH_UNITS = 'meters'
TIME_UNITS = 'days'

dirs = Dirs(HOME)

sim_name = 'MoervaarDpr'
dirs = dirs.add_case(sim_name)
os.chdir(dirs.case)
