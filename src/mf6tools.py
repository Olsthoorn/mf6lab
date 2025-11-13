# MFLAB file to read the parameters from Excel as used by mflab.py and mf_adapt.py

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class Dirs():
    """
    Class determining directory structure of the current model case
    
    General mf6 directory structure:
    mf6lab
        |__bin
        |__Projects
        |   |__project1
        |   |__project2
        |   |   |__cases
        |   |   |   |__case1
        |   |   |   |   |__data (case level data)
        |   |   |   |   |__doc (case level doc)
        |   |   |   |   |__GWF (modflow flow)
        |   |   |   |   |__GWT (modflow transport)
        |   |   |   |   |__images (case level images)
        |   |   |   |   |__MP7 (modpath)
        |   |   |   |   |__notebooks (case level notebooks)
        |   |   |   |   |__SIM (modflow sim)
        |   |   |   |   |__src (case level python source code)        
        |   |   |   |__case2
        |   |   |__data (project level data)
        |   |   |__doc (project level doc)
        |   |   |__images (project level images)
        |   |   |__notebooks (project level notebooks)
        |   |   |__src (project level source code)
        |   |   |__.code-workspace (project level workspace settings)
        |   |__project3
        |__src (mflab level source python code)
        .git (mf6 level .git)
        .venv (mf6 level virtual environment(s))
        .env (mf6 level .env (only read by debugger and interactive))
        .vscode (mf6 level settings)
           |__settings.json
           |__launch.json
        .gitignore (mf6 level .gitignore)
        LICENSE
        README.md

    
    HOME  # Home directory for this project
        .git   # Git directory should also link to github
        bin
            mf6.mac # Runnable code (hard linked to original)
            ..      # Other runnable code like Zonebudget
        cases   # folder collecting all cases for this projec
            <case_name1>  # First case <..> denote dedicated name.
                data
                doc
                GWF (= ground water flow)
                    mfsim.lst
                    ..
                    <case_name1>.lst
                    <case_name1>.dis
                    ..
                GWT (= ground water transport)
                     mfsim.lst
                    ..
                    <case_name1>.lst
                    <case_name1>.dis
                    ..
                images
                meteo
                MP7 (= Modpath)
                photos
                SIM (= )
                __init__.py
                mf_adapt.py
                mf_analyze.py
                <case_name1>.xlx
            
            <case_name2>  # Second case for this project
            <case_name3>  # Third case for this project
            ..
        data  # Data for this projec
            gis  # Shapefiles etc. (GIS) for this project
                <shape_file>.shp
        doc # Documentation for this project
            <documentation>.lyx
        images
        literature # General backgound references.
        lyx
        notebooks
        src     # Specific for this project
            __init__.py
            <fname1>.py
            <fname2>.py
            ..
    Individual cases use the <sim_name> and are placed under cases.
    Each case directory will contain the case-specific python source
    files among which are
        `mf_adapt.py`
        `mf_analyze.py`
        'settings.py' (if present imported by mf_adapt)
        and the parameter excel file f'{sim_name}.xlsx'
           
    @TO 231112, 25111
    """

    
    def __init__(self, HOME=None):
        """Return project directories class object.
        
        HOME must be the path to the project directory holding all the cases
        or None of VScode is launched from the project directory
        
        Parameters
        ----------
        home: str or path, default '.'
            home directory for this modelling projecta containing all cases.
        """
        if HOME is not None:
            HOME = Path(os.getcwd())
        MF6LAB = HOME.parent.parent

        if not str(MF6LAB).endswith('mf6lab'):
            raise ValueError(
                """`home` must the project directory
                or None if you launch VScode from the project directory.
                Not: {HOME}"""
                )
        
        self.HOME = HOME
        self.mf6lab = MF6LAB
        print(f"Project HOME directory: {HOME}")

        # --- Subdirectories
        self.bin = os.path.join(self.HOME, 'bin')
        self.cases =os.path.join(self.HOME, 'cases')
        self.data = os.path.join(self.HOME, 'data')
        self.gis = os.path.join(self.data, 'gis')
        self.doc = os.path.join(self.HOME, 'doc')
        self.images = os.path.join(self.doc, 'images')
        self.literature = os.path.join(self.doc, 'literature')
        self.src = os.path.join(self.HOME, 'src')
        
        # --- Make sure the src directory is known to Python
        sys.path.insert(0, self.src)
        
        # --- Create the directory if it doesn't yet exist
        self.create()
        
        
    def create(self):
        for d in [self.HOME, self.bin,
                  self.cases,
                  self.data,
                  self.gis,
                  self.doc,
                  self.images,                  
                  self.src]:
            if not os.path.isdir(d):
                os.mkdir(d)
            
    def add_case(self, case_name):
        """Add a case with subdirs and return the update dirs object.
        
        Parameters
        ----------
        case_name: str
            name of the case to be added to self.cases.
            If necessary the folder will be created.
            The path to the case is retured.
            self.f'{case_name}' can be used to address the folder.
        """
        self.case = os.path.join(self.cases, case_name)
        if not os.path.isdir(self.case):
            os.mkdir(self.case)
            print(f"dirs.case = '{self.case}' now exists")
            print(f"dirs.{case_name} = '{self.case}' now exists.")
        else:
            print(f"dirs.case = '{self.case}' already exists")
            print(f"dirs.{case_name} = '{self.case}' already exists.")
        assert os.path.isdir(self.case), f"Directory '{self.case}' does not exist."
            
        # --- Add subdirs to case directory if not already present
        for subd in ['SIM', 'GWF', 'GWT', 'MP7', 'data', 'doc', 'images', 'src']:
            subdir =  os.path.join(self.case, subd)           
            exec(f"self.{subd} = '{subdir}'")
            if not os.path.isdir(subdir):
                os.mkdir(subdir)                
                print(f"dirs.{subd} = '{subdir}' now exists.")
            else:
                print(f"dirs.{subd} = '{subdir}' already exists.")
            assert os.path.isdir(subdir), f"Directory '{subdir}' does not exist"
        return self
    
    def add_to_home(self, dirname):
        """Add a directory to the current HOME directory (self.HOME)."""
        folder = os.path.join(self.HOME, dirname)
        if os.path.exists(folder):
            print(f"{folder} already exists.")
        else:
            os.mkdir(folder)
            print(f"{folder} now exists.")
        exec(f"self.{dirname.replace('/','_')} = '{folder}'")
        
    def add_to_case(self, dirname):
        """Add a directory to the current project directory."""
        folder = os.path.join(self.case, dirname)
        if os.path.exists(folder):
            print(f"{folder} already exists.")
        else:
            os.mkdir(folder)
            print(f"{folder} now exists.")
        exec(f"self.{dirname} = '{folder}'")
    
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
    folder = '/Users/Theo/GRWMODELS/python/Pennink_model'
    f = os.path.join(folder, 'cases', 'Series2')
    os.chdir(f)

    sim_name = 'Series2'
    wbk_name = os.path.join(f, sim_name + '.xlsx')
    
    use_models, use_packages  = get_models_and_packages_from_excel(wbk_name, sheet_name='NAM')
    
    perioddata = get_periodata_from_excel(wbk_name, sheet_name='PER')
    
    print(sim_name)
    
    
    
    