# -*- coding: utf-8 -*-
"""_summary_

Settings (this file) can be run on its own, and will be imported
by mf_adapt to start running Modflow.

It specifies:
    properties of the simulation
    model properties not already in the origonal data (shpaefile)
"""

props = {
        'nper_test': 1000, # max number of stress periods when test is True
        'test': False,
        'length_units': 'meters',
        'time_units': 'days',
        'use_w_not_c': False,     # Use anal. form. instead of entry resistance and circumference.
        'start_date_time': '2024-01-27',
         'oc_frequency': 1,
         'icelltype': 1,
         'dx':   1.0,    # [m] cell width
         'minDz': 0.01, # m min layer thickness (also for pinched-out layers)
         'drain_depth': 0.15, # m
         'cDrainage':  100.0, # d
         'rch':        0.001, # m/d
         'strthd':       0.0, # m initial head
}

