# MFLAB file to read the parameters from Excel as used by mflab.py and mf_adapt.py

import os
import sys

import numpy as np
import pandas as pd
from pathlib import Path
import logging
from logging_setup import configure_logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
configure_logging()

class Dirs():
    """
    Class that generates the namespace for the current case
    
    General mflab directory tree:
    
    mf6lab

    |__bin
    |   |__mf6 executable (or symlink to it)
    |__doc (mf6level documentation)
    |__Projects
    |   |__<project1>
    |   |__<project2>
    |       |__cases
    |           |__<case1>
    |           |__<case2>
    |               |__data (case data)
    |                   |__GIS (case level GIS files)
    |                   |__meteo
    |                   |__<other data files>
    |               |__doc (case documentation)
    |               |__GWF (Modflow gwf files)
    |               |__GWT (Modflow transport files)
    |               |__images
    |               |__MP7 (Moddpath files)
    |               |__notebooks (jupyter .ipynb files)
    |               |__SIM (Modflow sim files)
    |               |__src (case-level .py files)
    |           |__<case3>
    |       |__data (project level data)
    |       |__doc (project level documentation)
    |       |__notebooks (project level notebooks)
    |       |__src (project level .py files)
    |   |__<project3>
    |__src (mf6lab level .py files)
    |__LICENSE
    |__README.MD

 
    @TO 20231112, 20251115
    """ 
    def __init__(self):
        """Return case directories class/namespace object creating directories if necessary.
        
        Launch dirs = Dirs from the case directory to get namespace of current case.

        Missing directories will be created.

        mf6lab/src, project/src and case/src will be added to sys.path
        
        >>>dirs = Dirs()
        >>>dirs.HOME
        >>>dirs.doc
        >>>dirs.mflab
        >>>dirs.case
        >>>dirs.src
        >>>os.path.join(dirs.case, 'src')
        >>>os.path.join(dirs.mf6lab, 'src')
        >>>for p in sys.path: print(p)                
        """
        ccd = Path(os.getcwd()) # Current Case Directory

        parts = ccd.parts
        try:
            idx = parts.index('mf6lab')
        except ValueError:
            raise FileNotFoundError(f"File not in mf6lab {ccd}")
        
        if len(parts) < idx + 5:
            raise ValueError(f"Script not lauchned from case sub folder {os.getcwd()}")

        self.mf6lab = str(os.path.join(*parts[:idx + 1], 'src'))
        self.proj   = str(os.path.join(*parts[:idx + 3], 'src'))
        self.case   = str(os.path.join(*parts[:idx + 5], 'src'))
        self.HOME   = self.case
        self.home   = self.case
        
        # --- make these three src directories known to python
        # --- already accomplished by mf6.bootstrap, so skip
        #for pth in [self.mf6lab, self.proj, self.case]:
        #    sys.path.insert(0, os.path.join(pth, 'src'))
            
        # --- Subdirectories for case
        case_dir = os.path.join(*Path(self.case).parts[:-1])
        subs = ['data', 'doc', 'GIS', 'GWF', 'GWT', 'meteo', 'MP7', 'images', 'notebooks', 'SIM', 'src']        
        for sub in subs:
            pth = os.path.join(case_dir, sub)
            setattr(self, sub, pth)
            if not os.path.isdir(pth):
                os.mkdir(pth)
                print(f"folder {pth} created")
            else:
                print(f"folder {pth} already exists")
    
    def __str__(self):
        return self.__doc__


def get_models_and_packages_from_excel(wbk_name, sheet_name='NAM'):
    """Return which models and packages will be used in the simulation.
  
    Parameters
    ----------
    wbk_name: str
        Excel workbook name that holds the mf6 parameters in sheet 'MF6'
    sheet_name: str
        The case-sensitive name of the sheet in the workbook telling
        which models and packages will be used. Should always be 'NAM'
        
    Returns
    -------
    dict with models and packages to be used.
    """
    # --- Packages is from line with 'Package'
    packages = pd.read_excel(wbk_name, sheet_name=sheet_name, header=0,
                                index_col='ModelPkg')

    # --- Only need the first column when ON / OFF is True
    packages = list(packages.index[packages['ON / OFF'] > 0])
    
    # --- Return all lower case with first letter capitalized
    packages = [m[0:1].upper() + m[1:].lower() for m in packages]
    models = list(np.unique([p[:3] for p in packages])) # keeps order i n list
    return models, packages
    
class ExecValue:
    """Class to convert quoted string to list or tuple (inloop exec)."""
    def __init__(self, value):
        exec(f"self.value = {value}")


def get_mf6_params_from_excel(wbk_name, sheet_name='GWF6'):
    """Read mf6 parameters from workbook (sheetname='GWF6').
    
    Read the parameters from Excel workbook as a pd.DataFrame and turn
    them into a dictionary, where the keys are Modlow package names.
    
    The headers in the sheet of the workbook are:
    ['Package', 'Param', 'Value']
    
    Parameters
    ----------
    wbk_name: str
        Excel workbook name that holds the mf6 parameters in sheet 'MF6'
    sheet_name: str
        The case-sensitive name of the sheet in the workbook.
        Examples are MF6, PER and LAY
        
    Returns
    -------
    params: dictionary or pd.DataFrame
        a dict with package name as key, in which each item is a dictionary
        specifying the values for this key.
        
    @TO 220413
    """
    paramsDf = pd.read_excel(wbk_name, sheet_name=sheet_name,
                             header=0, usecols="A:C", engine="openpyxl"
                             ).dropna(axis=0)
    
    params = dict()
    
    # --- Each package gets its own subdict with the actual values
    for pkg in np.unique(paramsDf[paramsDf.columns[0]].values):
        params[pkg]=dict()
        
    # --- Convert package parameters to correct type using the column "type"
    for i in paramsDf.index:
        pkg, param, value = paramsDf.loc[i]
        if isinstance(value, str):
            if value  == 'None':
                value = None
            elif value[0] in '[()':
                value = ExecValue(value).value
            elif "'" in value:
                value = value.replace("'","")
            elif '"' in value:
                value = value.replace('"', '')
            else:
                pass
        elif isinstance(value, (bool, float)):
            pass
        elif isinstance(value, int):
            if param[0] in 'iIjJkKlLmMnN':
                value = float(value)
            else:
                pass
        else: # --- Immediately verify unknown types
            raise ValueError("Unknown parameter type pkg={}, {}, {}, {}".format(
                            pkg, param, value, type))
        params[pkg][param]=value        
    return params


def get_periodata_from_excel(wbk_name, sheet_name='PER'):
    """Return a pd.DataFrame with the periodata read from excel workbook.
    
    Parameters
    ----------
    wbk_name: str
        name of the parameter workbook (Excel file)
    sheet_name: str, default 'PER'
        the (case senitive) name of the worksheet, 'PER'
        Notice that IPER in the worksheet must be zero-based.
        
    Returns
    -------
    period_data
        pd.DataFrame with index IPER and columns PERLEN, NSTP, TSMULT and others.
        The index is filled up in front of the numbered IPER.
        So if the first IPER is 7 then lines 0-7 will be filled with the values
        of this line. If the next IPER = 12, then lines 8 - 12 will be filed with the
        line with IPER = 8, etc. Hence, the last line is the last period number and
        the total number of stress periods is one more than this.
    """
    p_data = pd.read_excel(wbk_name, sheet_name=sheet_name, header=1,
                                index_col='IPER')
    p_data.index = p_data.index.astype(int)
    
    # --- Remove any dummy stress period lines (index <= 0)
    p_data = p_data.loc[p_data.index >= 0]

    # Set type of stress period data
    p_data = p_data.astype({'PERLEN': float, 'NSTP': int, 'TSMULT': float})
    
    p_data = p_data.loc[p_data['PERLEN'] > 0]
    
    # --- IPER (=index) in PER sheet are zero based, keep it that way.
    
    p_index = np.asarray(p_data.index)
    period_data = pd.DataFrame(index=np.arange(p_index[-1] + 1), columns=p_data.columns)
    
    # Forward fill period_data
    ip = 0
    for idx in p_index:
        while ip <= idx:
            period_data.loc[ip] = p_data.loc[idx]
            ip += 1
   
    period_data['NSTP'] = period_data['NSTP'].astype(int)
    return period_data

# Same for PER and LAY


def show_animation_progress(frame, nbreak=50, ntot=None):
    """"Show progress of animation.
    
    Parameters
    ----------
    frame: int
        animation frame
    nbreak: int
        Print "frame counts and newline each nbreak frames.
    ntot: int
        total number of frames in the animation.
        Print frame count and newline after last frame was done.
    """
    print('.', end="")
    if (frame + 1) % nbreak == 0 or (frame + 1) == ntot:
        print('{} frames\n'.format(frame + 1))


if __name__ == '__main__':
  
    ccd = Path('/Users/Theo/GRWMODELS/python/mf6lab/Projects/Pennink-Model/cases/Series1').resolve()
    os.chdir(ccd)
    
    print()
    print('ccd (current case directory) : ', ccd)
    
    dirs = Dirs()

    print("-------------")
    print(dirs.mf6lab)
    print(dirs.proj)
    print(dirs.case)
    print(dirs.data)
    print(dirs.home)
    print(dirs.images)
    print(dirs.doc)
    print("-------------")

    case_name = Path(dirs.case).parts[-1]
    wbk_name = os.path.join(dirs.case, case_name + '.xlsx')
    
    assert os.path.isfile(wbk_name), f"No file {wbk_name}"
    
    use_models, use_packages  = get_models_and_packages_from_excel(wbk_name, sheet_name='NAM')
    
    perioddata = get_periodata_from_excel(wbk_name, sheet_name='PER')
    
    print("Case_name: ", case_name)
    print("wbk_name : ", wbk_name)
    print()
     
    
    
    