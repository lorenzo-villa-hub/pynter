
import numpy as np
from scipy.optimize import bisect, root_scalar, minimize_scalar, minimize
from monty.json import MontyDecoder, MSONable, MontyEncoder
from abc import ABCMeta
import pandas as pd
import os
import os.path as op
import json
import copy
import warnings
import matplotlib.pyplot as plt

from pymatgen.electronic_structure.dos import FermiDos

from .chempots.oxygen import get_pressure_reservoirs_from_precursors, get_oxygen_pressure_reservoirs
from .corrections.kumagai import get_kumagai_correction
from .corrections.freysoldt import get_freysoldt_correction_from_locpot
from .defects import Defect, get_defect_from_string
from .electronic_structure import get_carrier_concentrations
from .entries import DefectEntry, _get_computed_entry_from_path
from .plotter import (
                    plot_pO2_vs_concentrations,
                    plot_variable_species_vs_concentrations,
                    plot_binding_energies,
                    plot_formation_energies,
                    plot_charge_transition_levels
                    )
from pynter.tools.utils import get_object_feature, select_objects, sort_objects


class DefectsAnalysis(MSONable,metaclass=ABCMeta):
    
    def __init__(self, entries, vbm, band_gap, sort_entries=True):
        """ 
        Class to compute collective properties starting from single calculations of point defects.

        Parameters
        ----------
        entries: (list) 
            A list of DefectEntry objects.
        vbm: (float) 
            Valence band maximum of the pristine material in eV.
        band_gap : (float)
            Band gap of the pristine material in eV
        """
        self.entries = self.sort_entries(inplace=False,entries=entries) if sort_entries else entries
        self.vbm = vbm
        self.band_gap = band_gap
        self.groups = self._group_entries()
        self.names = list(self.groups.keys())
        self._thermodata = None
        

    def __str__(self):     
        return self.table().__str__()
            
    def __repr__(self):
        return self.table().__repr__()
    
    def _repr_html_(self):
        return self.table()._repr_html_()
    
    def __iter__(self):
        return self.entries.__iter__()
    
    def __getitem__(self,index):
        return self.entries.__getitem__(index)
    
    def append(self,item):
        return self.entries.append(item)
    
    def _group_entries(self):
        groups = {}
        for e in self.entries:
            if e.name not in groups.keys():
                groups[e.name] = [e]
            else:
                groups[e.name].append(e)
        return groups
            
    def copy(self):
        return DefectsAnalysis(self.entries, self.vbm, self.band_gap)
    
    @property
    def thermodata(self):
        """
        Result of thermodynamics calculation (`plot_brouwer_diagram` or `plot_doping_diagram functions`).
        """    
        return self._thermodata
    
    @property
    def elements(self):
        """
        List of all the elements involved in exchange of atoms with a reservoir 
        (elements present in entry.delta_atoms for all entries).
        """
        elements = []
        for entry in self.entries:
            for el in entry.delta_atoms.keys():
                if el not in elements:
                    elements.append(el)
        return elements


    def as_dict(self):
        """
        Returns:
            Json-serializable dict representation of DefectsAnalysis
        """
        d = {
        "entries" : MontyEncoder().encode(self.entries),
        "vbm":self.vbm,
        "band_gap":self.band_gap
            }
        return d


    @classmethod
    def from_dict(cls,d):
        """
        Reconstruct a DefectsAnalysis object from a dict representation created using
        as_dict().

        Parameters
        ----------
            d (dict): dict representation of DefectsAnalysis.

        Returns:
        -------
            DefectsAnalysis object
        """
        entries = MontyDecoder().decode(d['entries'])        
        vbm = d['vbm']
        band_gap = d['band_gap']
        return cls(entries,vbm,band_gap)


    @staticmethod
    def from_dataframe(df,vbm, band_gap):
        """
        Create DefectsAnalysis object from pandas DataFrame (df). If df has been 
        exported with include_structure=True, DefectEntry objects will be 
        imported including defect and bulk structures.

        Parameters
        ----------
        df : (DataFrame)
            pandas DataFrame. Needs to be in the following format:

            - rows are defect entries

            - non-optional columns are:
                - name : defect name in the format given by Defect objects. 
                - charge : defect charge
                - multiplicity : site multiplicity of the defect
                - energy_diff : difference in energy btw defect and bulk structures

            - If the df was not exported from this class, you need to provide also:
                - bulk_volume : Volume of the bulk cell in A°3(unit cell or supercell 
                                relative to your multiplicity)
        vbm : (float)
            Valence band maximum of bulk structure in eV.
        band_gap : (float)
            Band gap of bulk structure in eV.
        
        Returns
        -------
        DefectsAnalysis object
        """
        default_columns = ['name','charge','multiplicity','energy_diff','bulk_volume']
        entries = []
        for idx,row in df.iterrows():
            if 'defect' in df.columns:
                defect = row['defect']
            else:
                defect = get_defect_from_string(
                                        row['name'],
                                        charge=row['charge'],
                                        multiplicity=row['multiplicity'],
                                        bulk_volume=row['bulk_volume']
                                        )
            energy_diff = row['energy_diff']
            corrections = {}
            data = {}
            for col in df.columns:
                if 'corr' in col:
                    key = col.split('_')[1]
                    corrections[key] = row[col]
                elif col not in default_columns:
                    if col not in DefectEntry.__dict__.keys(): # check data is not already a property
                        data[col] = row[col]

            entry = DefectEntry(
                            defect=defect,
                            energy_diff=energy_diff,
                            corrections=corrections,
                            data=data)
            entries.append(entry)
        
        return DefectsAnalysis(entries=entries,vbm=vbm,band_gap=band_gap)
    

    @staticmethod
    def from_file(filename,vbm,band_gap,format=None, **kwargs): 
        """
        Create DefectsAnalysis object from file.
        Available formats are:
            - 'json'
            - 'pkl'
            - 'csv'
        Check docs in `from_dataframe` function for file formatting requirements.
        """
        if '.json' in filename or format == 'json':
            return DefectsAnalysis.from_json(filename)
        elif '.pkl' in filename or format == 'pkl':
            df = pd.read_pickle(filename, **kwargs)
        elif '.csv' in filename or format == 'csv':
            if kwargs:
                if 'index_col' not in kwargs:
                    kwargs['index_col'] = False
            else:
                kwargs['index_col'] = False
            df = pd.read_csv(filename, **kwargs)
        else:
            raise ValueError('Invalid file format, available are "json","pkl","csv" and "excel"')
        return DefectsAnalysis.from_dataframe(df=df,vbm=vbm,band_gap=band_gap)


    @staticmethod
    def from_json(path_or_string):
        """
        Build DefectsAnalysis object from json file or string.

        Parameters
        ----------
        path_or_string : (str)
            If an existing path to a file is given the object is constructed reading the json file.
            Otherwise it will be read as a string.

        Returns
        -------
        DefectsAnalysis object.
        """
        if op.isfile(path_or_string):
            with open(path_or_string) as file:
                d = json.load(file)
        else:
            d = json.loads(path_or_string)
        return DefectsAnalysis.from_dict(d)


    @staticmethod
    def from_vasp_directories(
                            path_defects,
                            path_bulk,
                            common_path=None,
                            get_charge_correction='kumagai',
                            dielectric_tensor=None,
                            get_multiplicity=False,
                            get_data=True,
                            band_gap=None,
                            vbm=None,
                            initial_structure=False,
                            function=None,
                            computed_entry_kwargs={},
                            finder_kwargs={},
                            correction_kwargs={}):
        """
        Generate DefectsAnalysis object from VASP directories read with Pymatgen.

        Parameters
        ----------
        path_defects : (str)
            Path of VASP defects calculation.
        path_bulk : (str)
            Path of VASP bulk calculation. 
        get_charge_correction : (str or bool)
            Compute charge corrections from VASP directories. 
            To skip corrections set it to False.
            Methods available are:
                - "kumagai" : Extended FNV scheme (eFNV).
                - "freysoldt" : FNV scheme.
        get_multiplicity : (bool)

        get_data : (True)
            Store path in DefectEntry.data dictionary.
        band_gap : (float)
            Band gap of bulk material. If not provided is determined from bulk VASP directory.
        vbm : (float)
            Valence band maximum of bulk material. If not provided is determined from bulk VASP directory.
        initial_structure : (bool)
            Use initial structure for defect recognition. Useful when relaxations are large and 
            defect_finder struggles to find the right defects.
        function : (func)
            Function to apply to each DefectEntry. Useful to automate custom entry modification.
        computed_entry_kwargs : (dict)
            Kwargs to pass to `Vasprun.get_computed_entry`.
        finder_kwargs : (dict)
            Kwargs to pass to `defect_finder`.
        correction_kwargs : (dict)
            Kwargs to pass to the charge correction methods.

        Returns
        -------
        DefectsAnalysis object.
        """ 
        from pymatgen.io.vasp.outputs import Vasprun
        if band_gap and vbm:
            parse_eigen = False
        else:
            parse_eigen = True
        vasprun_bulk = Vasprun(op.join(path_bulk,'vasprun.xml'),parse_dos=False,parse_eigen=parse_eigen,parse_potcar_file=False)
        computed_entry_bulk = vasprun_bulk.get_computed_entry(**computed_entry_kwargs)
        band_gap = band_gap or vasprun_bulk.eigenvalue_band_properties[0]
        vbm = vbm or vasprun_bulk.eigenvalue_band_properties[2]

        common_path = common_path or path_defects # if not set enter all VASP directories
        entries = []
        for root,dirs,files in os.walk(path_defects):   
            if 'vasprun.xml' in files and common_path in root:
                path = op.abspath(root)
                print('importing from %s'%path)                

                data = {'path':path} if get_data else None
                entry = DefectEntry.from_vasp_directories(
                                                        path_defect=path,
                                                        computed_entry_bulk=computed_entry_bulk,
                                                        corrections={}, # get corrections and multiplicity automatically?
                                                        multiplicity=1,
                                                        data=data,
                                                        label=None,
                                                        function=function,
                                                        initial_structure=initial_structure,
                                                        computed_entry_kwargs=computed_entry_kwargs,
                                                        finder_kwargs=finder_kwargs)
                entries.append(entry)

                ck = correction_kwargs
                ck['charge'] = entry.charge
                if get_charge_correction == 'kumagai':
                    ck['defect_path'] = path
                    ck['bulk_path'] = path_bulk
                    if dielectric_tensor:
                        ck['dielectric_tensor'] = dielectric_tensor
                    ck['get_correction_data'] = False
                    corr = get_kumagai_correction(**ck)

                elif get_charge_correction == 'freysoldt':
                    defect_path_locpot = op.join(path,'LOCPOT')
                    bulk_path_locpot = op.join(path_bulk,'LOCPOT')
                    ck['defect_path_locpot'] = defect_path_locpot
                    ck['bulk_path_locpot'] = bulk_path_locpot
                    ck['finder_kwargs'] = finder_kwargs
                    if dielectric_tensor:
                        ck['dielectric_tensor'] = dielectric_tensor
                    corr = get_freysoldt_correction_from_locpot(**ck)

                if get_charge_correction:
                    if type(corr) == tuple:
                        corr = corr[0]
                        plt.show()
                    entry.set_corrections(**{get_charge_correction:corr})
        
        return DefectsAnalysis(entries=entries,vbm=vbm,band_gap=band_gap)



    def to_dataframe(self,
                    entries=None,
                    include_structures=False,
                    include_data=True,
                    properties=[],
                    functions={}):
        """
        Export DefectsAnalysis object as DataFrame. 
        Default exported columns are:
            - "name"
            - "charge"
            - "multiplicity"
            - "energy_diff"
            - "bulk_volume"
        
        Corrections dict (`entry.corrections`) are exported as separate columns
        for each keys in this format "corr_{key}".

        Parameters
        ----------
        entries : (list)
            Entries to export. If None all entries are displayed.
        include_structures : (bool)
            Add full Defect object to DataFrame containing defect and bulk structures.
            This makes the DataFrame less portable. 
        include_data: (bool)
            Include entry.data dictionary. 
        properties: (list)
            List of strings with defect entry attributes or methods to include in df.
        functions : (dict)
            Dictionary with column names as keys and functions as values. 
            The function needs to take a DefectEntry as input and the output
            will be stored in df.
        Returns
        -------
        df : 
            DataFrame object.
        """
        if not entries:
            entries = self.entries
        d = {}
        index = []
        table = []
        for e in entries:
            symbol = e.symbol
            d = {}
            d['name'] = e.name
            d['charge'] = e.charge
            d['multiplicity'] = e.multiplicity
            d['energy_diff'] = round(e.energy_diff,4)
            if include_structures:
                d['defect'] = e.defect
            if e.corrections:
                for key,value in e.corrections.items():
                    d[f'corr_{key}'] = value
            d['bulk_volume'] = round(e.defect.bulk_volume, 4)
            if include_data:
                if e.data:
                    for key,value in e.data.items():
                        d[key] = value
            if properties:
                for feature in properties:
                    if isinstance(feature,list):
                        raise ValueError('Only simple attributes can be exported')
                    else:
                        key = feature
                    d[key] = get_object_feature(e,feature)
            if functions:
                for key, fn in functions.items():
                    d[key] = fn(e)        

            table.append(d)
        df = pd.DataFrame(table)

        return df
    

    def to_file(self,filename,format=None,export_kwargs={},**kwargs):
        """
        Export DefectsAnalysis object to file. 

        Parameters
        ----------
        filename : (str)
            Path of file. 
        format : (str)
            If not specified the file extension in `filename` is used.
            Formats available:
            - "pkl" : pickle file. Allows to store structures (default).
            - "json" : DefectsAnalysis object as json file.
            - "csv" : Does not allow to store structures. 
        export_kwargs : (dict)
            Kwargs to pass to file exporting function
        kwargs : (dict)
            Kwargs to pass to to_dataframe function. If not provided, 
            'include_strctures' is to True if chosen format is 'pkl', 
            set to False if chosen format is 'csv' or 'excel'.
        """
        if format == 'json' or '.json' in filename:
            self.to_json(path=filename) 
        elif format == 'pkl'  or '.pkl' in filename:
            if 'include_structures' not in kwargs.keys():
                kwargs['include_structures'] = True
            df = self.to_dataframe(**kwargs)
            df.to_pickle(filename,**export_kwargs)
        elif format == 'csv' or '.csv' in filename:
            if 'include_structures' not in kwargs.keys():
                kwargs['include_structures'] = False
            df = self.to_dataframe(**kwargs)
            if export_kwargs:
                if 'index' not in export_kwargs:
                    export_kwargs['index'] = False
            else:
                export_kwargs = {'index': False}
            df.to_csv(filename,**export_kwargs)


    def to_json(self,path='',**kwargs):
        """
        Save DefectsAnalysis object as json string or file.
        Parameters
        ----------
        path : (str), optional
            Path to the destination file. If '' the path is set to "self.path/self.name.json".
            If None a string is exported. 
        Returns
        -------
        d : (str)
            If path is not set a string is returned.
        """
        self.reset_all_custom_functions()
        d = self.as_dict(**kwargs)
        if path == '':
            path = op.join(self.path,self.name+'.json')
        if path:
            with open(path,'w') as file:
                json.dump(d,file)
            return
        else:
            return d.__str__() 
        



    def binding_energy(self,name,fermi_level=0, temperature=0, **eform_kwargs):
        """
        Compute the binding energy for a defect complex as:
        Eb = Ef(complex) - \sum_D Ef(D), where D are the individual defects.
        
        Parameters
        ----------
        name :(string)
            Name of defect complex as assigned in defect entry object
        fermi_level : (float)
            Fermi level in eV.
        temperature : (float)
            Temperature in K. If no custom formation energy is provided, this arg has no effect.
        eform_kwargs : (dict)
            Kwargs to pass to `entry.formation_energy`.

        Returns:
        --------
            float : Binding energy in eV
        """
        stable_charges = self.stable_charges(
                                    chemical_potentials=None,
                                    fermi_level=fermi_level,
                                    temperature=temperature,
                                    entries=None,
                                    **eform_kwargs)
        
        charge = stable_charges[name][0]
        binding_energy = stable_charges[name][1]
        entry = self.select_entries(names=[name],charge=charge)[0]
        for sd in entry.defect.defects: #subtract energies of single defects
            dname = sd.name
            binding_energy = binding_energy - stable_charges[dname][1]    
        
        return binding_energy


    def carrier_concentrations(self,bulk_dos,temperature=300,fermi_level=0.):
        """
        Get intrinsic carrier concentrations by integrating over the density of states
        given a fixed Fermi level.

        Parameters
        ----------
            dos : (dict or Dos)
                Density of states to integrate. Can be provided as density of states D(E)
                or using effective masses.
                Format for effective masses:
                    dict with following keys:
                        - "m_eff_h" : holes effective mass in units of m_e (electron mass)
                        - "m_eff_e" : electrons effective mass in units of m_h          
                        - `band_gap` needs to be provided in args
                Formats for explicit DOS:
                    dictionary with following keys:
                        - 'energies' : list or np.array with energy values
                        - 'densitites' : list or np.array with total density values
                        - 'structure' : pymatgen Structure of the material, 
                                        needed for DOS volume and charge normalization.
                    or a pymatgen Dos object (Dos, CompleteDos or FermiDos).

            fermi_level : (float)
                The Fermi level relative to the VBM in eV.
            temperature : (float)
                The temperature in Kelvin.

        Returns:
        --------
            h : (float)
                Absolute value of hole concentration in 1/cm^3
            n : (float)
                Absolute value of electron concentration in 1/cm^3
        """   
        return get_carrier_concentrations(
                                    dos=bulk_dos,
                                    fermi_level=fermi_level,
                                    temperature=temperature,
                                    band_gap=self.band_gap)


    def charge_transition_levels(self,
                                energy_range=None,
                                temperature=0,
                                entries=None,
                                get_integers=True,
                                **eform_kwargs):
        """
        Computes charge transition levels for all defect entries.
        
        Parameters
        ----------
        energy_range : (list or tuple)
            Energy range in eV to evaluate the charge transition levels, default to (-0.5, band_gap + 0.5).   
        temperature : (float)
            Temperature in K. If no custom formation energy is provided, this arg has no effect.
        entries : (list)
            List of entries to calculate. If None all entries are considered.
        get_integers : (bool)
            Save charges as integers. More convenient for plotting.
        eform_kwargs : (dict)
            Kwargs to pass to `entry.formation_energy`.
        
        Returns:
        --------
        Dictionary with defect name and list of tuples for charge transition levels in eV:
                {name:[(q1,q2,ctl),(q2,q3,ctl),...}
        """                
        entries = entries if entries else self.entries
        if energy_range == None:
            energy_range = (-0.5,self.band_gap +0.5) 
        # creating energy array
        npoints = 3000
        step = abs(energy_range[1]-energy_range[0])/npoints
        e = np.arange(energy_range[0],energy_range[1],step)
        
        # starting point is first E value
        previous_stable_charges = self.stable_charges(
                                    chemical_potentials=None,
                                    fermi_level=energy_range[0],
                                    temperature=temperature,
                                    entries=entries,
                                    **eform_kwargs)
        
        charge_transition_levels = {name:[] for name in previous_stable_charges}
        for i in range(0,len(e)):
            stable_charges = self.stable_charges(
                                    chemical_potentials=None,
                                    fermi_level=e[i],
                                    temperature=temperature,
                                    entries=entries,
                                    **eform_kwargs)
            
            for name in stable_charges:
                new_charge = stable_charges[name][0]
                previous_charge = previous_stable_charges[name][0]
                if new_charge != previous_charge:
                    if name not in charge_transition_levels.keys():
                        charge_transition_levels[name] = []
                    previous_stable_charges = stable_charges
                    if get_integers:
                        previous_charge = int(previous_charge)
                        new_charge = int(new_charge)
                    charge_transition_levels[name].append((previous_charge,new_charge,e[i]))
                            
        return charge_transition_levels        
    

    def defect_concentrations(self,
                            chemical_potentials,
                            temperature=300,
                            fermi_level=0.,
                            fixed_concentrations=None,
                            per_unit_volume=True,
                            eform_kwargs={},
                            **kwargs):
        """
        Compute concentrations of all defects.

        For each entry, is possible to set a custom function to compute the defect 
        concentration using `set_defect_concentration_function`. If that is set,
        the new function is used instead of the default. 
        The function args must be the following:
            vbm : (float)
                Valence band maximum of bulk calculation in eV
            chemical_potentials : (dict)
                Chemical potentials of the elements involved in the defect.
            temperature : (float)
                Temperature in Kelvin.
            fermi_level : (float)
                Fermi level in eV (with respect to the VBM)
            per_unit_volume : (bool)
                Compute concentrations per unit volume using `self.defect.bulk_volume`.
            eform_kwargs : (dict)
                Kwargs to pass to `self.formation_energy`.
            kwargs : (dict)
                Additional custom kwargs.
        It is possible to set custom functions collectively using the 
        `set_defect_concentration_functions` method in this class.
        To reset functions to default, use `reset_defect_concentration_functions`.

        If fixed_concentrations is provided the concentration of defect entries are 
        corrected according to the fixed provided values. More details can be found in 
        https://doi.org/10.1103/PhysRevB.106.134101 .
            
        Parameters
        ----------
        chemical_potentials: (dict)
            Dictionary of chemical potentials ({element: chempot})   
        temperature: (float) 
            Temperature in K.
        fermi_level: (float) 
            Fermi level in eV relative to valence band maximum.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys can be simple element strings
            (or vacancies of elements in the format 'Vac_{el}') if that element needs to be 
            fixed across all defect species, alternatively defect entry names can be used as well 
            to target specific defect entries. The values are the concentrations.
        per_unit_volume: (bool)
            Get concentrations in cm^-3. If False they are per unit cell.
        eform_kwargs : (dict)
                Kwargs to pass to `self.formation_energy`.
        kwargs : (dict)
            Additional custom kwargs to pass to `entry.defect_concentration`.

        Returns:
        --------
        DefectConcentrations object, behaves like a list.
        """
        concentrations = []
        if fixed_concentrations:
            dc = self.defect_concentrations(
                                chemical_potentials=chemical_potentials,
                                temperature=temperature,
                                fermi_level=fermi_level,
                                fixed_concentrations=None,
                                per_unit_volume=per_unit_volume,
                                eform_kwargs=eform_kwargs,
                                **kwargs)
            
            frozen = fixed_concentrations 

        for e in self.entries:
            # frozen defects approach
            if fixed_concentrations:
                c = e.defect_concentration(
                                vbm=self.vbm,
                                chemical_potentials=chemical_potentials,
                                temperature=temperature,
                                fermi_level=fermi_level,
                                per_unit_volume=per_unit_volume,
                                eform_kwargs=eform_kwargs,
                                **kwargs)
                
                corr = self._get_frozen_correction(e,frozen,dc)
                c = c * corr
                defconc = SingleDefConc(name=e.name,charge=e.charge,conc=c)
                concentrations.append(defconc)     
            
            else:
                c = e.defect_concentration(
                                vbm=self.vbm,
                                chemical_potentials=chemical_potentials,
                                temperature=temperature,
                                fermi_level=fermi_level,
                                per_unit_volume=per_unit_volume,
                                eform_kwargs=eform_kwargs,
                                **kwargs)
                
                defconc = SingleDefConc(name=e.name,charge=e.charge,conc=c)
                concentrations.append(defconc)
            

        return DefectConcentrations(concentrations)


    def _get_frozen_correction(self,e,frozen,dc):
        corr = 1
        df = Defect.from_string(e.name)
        for defect in df:
            typ, specie, name = defect.type, defect.specie, defect.name
            if typ == 'Vacancy':
                k = name
                if k in frozen.keys():
                    eltot = dc.get_element_total(specie,vacancy=True)
                    corr = corr * (frozen[k]/eltot)
            else:
                if name in frozen.keys():
                    dtot = 0
                    for n,c in dc.total.items(): # sum over all defects containing the specific specie
                        if name in n:
                            dtot += c
                    corr = corr * (frozen[name]/dtot) 
                else:
                    k = specie
                    if k in frozen.keys():
                        eltot = dc.get_element_total(specie,vacancy=False)
                        corr = corr * (frozen[k]/eltot)
        return corr

    
    def filter_entries(self,
                       inplace=False,
                       entries=None,
                       mode='and',
                       exclude=False,
                       types=None,
                       elements=None,
                       names=None,
                       function=None,
                       **kwargs):
        """
        Filter entries based on different criteria. Return another DefectsAnalysis object.

        Parameters
        ----------
        inplace : (bool)
            Apply changes to current DefectsAnalysis object.
        entries : (list)
            List of defect entries.
        mode : (str), optional
            Filtering mode, possibilities are: 'and' and 'or'. The default is 'and'. 
        exclude : (bool), optional
            Exclude the entries satisfying the criteria instead of selecting them. The default is False.
        types : (list), optional
            Class name of the defect in the entry.
        elements : (list), optional
            List of symbols of elements that need to belong to the defect specie.
            If None this criterion is ignored. The default is None.
        names : (list)
            List of entry names.
        function : (function), optional
            Specific funtion for more complex criteria. The function must take a DefectEntry
            object as argument and return a bool.
        **kwargs : (dict)
            Properties that the jobs need to satisfy. Keys are referred to attributes/methods 
            of the defect entry. To address more than one condition relative to
            the same attribute, use lists or tuples (e.g. charge=[0,1]).

        Returns
        -------
        DefectsAnalysis object.
        """
        output_entries = self.select_entries(entries=entries,mode=mode,exclude=exclude,types=types,
                                             elements=elements,names=names,function=function,**kwargs)
        
        if inplace:
            self.entries = output_entries
            return
        else:
            return DefectsAnalysis(output_entries, self.vbm, self.band_gap)


    def formation_energies(self,
                           chemical_potentials=None,
                           fermi_level=0,
                           temperature=0,
                           entries=None,
                           **eform_kwargs):
        """
        Compute formation energies for all defect entries.
        Returns a dictionary with names of defect entries as keys and 
        a list of tuples (charge,formation_energy) as values.    

        For each entry, is possible to set a custom function to compute the formation 
        energy using `set_formation_energy_function`. If that is set,
        the new function is used instead of the default. 
        The function args must be the following:
            vbm : (float)
                Valence band maximum of bulk calculation in eV
            chemical_potentials : (dict)
                Chemical potentials of the elements involved in the defect.
            fermi_level : (float)
                Fermi level in eV (with respect to the VBM).
            temperature : (float)
                Temperature in Kelvin.
            kwargs : (dict)
                Additional custom kwargs.
        It is possible to set custom functions collectively using the 
        `set_formation_energy_functions` method in this class.
        To reset functions to default, use `reset_formation_energy_functions`.

        
        Parameters
        ----------
        chemical_potentials: (dict)
            Dictionary of chemical potentials ({element: chempot})   
        fermi_level: (float) 
            Fermi level in eV relative to valence band maximum.
        temperature : (float)
            Temperature in K. If no custom formation energy is provided, this arg has no effect.
        entries : (list)
            List of defect entries to calculate. If None all entries are considered.
        eform_kwargs : (dict)
            Additional custom kwargs to pass to `entry.formation_energy`

        Returns
        -------
        Dictionary in the format:
            {name: [(charge,formation energy)] }
        """
        formation_energies = {}
        entries = entries if entries else self.entries
        for entry in entries:
            name = entry.name
            charge = entry.charge
            eform = entry.formation_energy(
                                vbm=self.vbm,
                                chemical_potentials=chemical_potentials,
                                fermi_level=fermi_level,
                                temperature=temperature,
                                **eform_kwargs)
            
            if name in formation_energies:
                formation_energies[name].append((charge,eform))
            else:
                formation_energies[name] = []
                formation_energies[name].append((charge,eform))
        
        return formation_energies
        
           
    def get_charge_transition_level(self,name,q1,q2,temperature=0,**eform_kwargs):
        """
        Compute charge transition level.

        Parameters
        ----------
        name : (str)
            Defect name in the format given by the Defect object.
        q1 : (int or float)
            Charge state (in units of electron charge).
        q2 : (int or float)
            Charge state (in units of electron charge).
        temperature : (float)
            Temperature in K. If no custom formation energy is provided, this arg has no effect.
        eform_kwargs : (dict)
            Kwargs to pass to `entry.formation_energy

        Returns
        -------
        Charge transition level (float)
        """
        chemical_potentials = None
        formation_energies = self.formation_energies(
                                chemical_potentials=chemical_potentials,
                                fermi_level=0,
                                temperature=temperature,
                                **eform_kwargs)
        
        for d in formation_energies[name]:
            if d[0] == q1:
                E1 = d[1]
            if d[0] == q2:
                E2 = d[1]
        charge_transition_level = -1*(E1 - E2)/(q1 - q2)
        
        return charge_transition_level
        

    def merge_entries(self,*args):
        """
        Create a new DefectsAnalysis object merging the defect entries of other objects.
        The VBM and band gap of all objects need to match for consistency.
        """
        new_entries = self.entries
        for da in args:
            current_index = args.index(da)
            if current_index != 0:
                previous_da = args[current_index-1]
                if da.vbm != previous_da.vbm and da.band_gap != previous_da.band_gap:
                    raise ValueError('Cannot merge entries, band gap and VBM are different')
            new_entries = new_entries + da.entries
        return DefectsAnalysis(new_entries, self.vbm, self.band_gap)
        

    def plot_brouwer_diagram(self,
                            bulk_dos,
                            temperature,
                            quench_temperature=None,
                            quenched_species=None,
                            quench_elements=False,
                            fixed_concentrations=None,
                            external_defects=[],
                            reservoirs = None,
                            precursors = None,
                            oxygen_ref = None,
                            pressure_range=(1e-20,1e10),
                            npoints = 50,
                            xtol=1e-10,
                            eform_kwargs={},
                            dconc_kwargs={},
                            **kwargs):
        """
        Plot Brouwer diagram (defect concentrations vs oxygen partial pressure). Wrapper function for 
        `pynter.defects.thermodynamics.DefectThermodynamics` and `pynter.defects.plotter`.
        If you need more control use the classes individually. 

        If `quenched_temperature` is set, defect concentrations are computed at the initial `temperature`, 
        and defect quilibrium in solved at `quenched_temperature`. For the same defect species, the charge
        state follow the temperature and Fermi level dependency at `quench_temperature`, the total
        concentration is kept fix at the initial `temperature`.
        If additional `fixed_concentrations` are provided, the concentration of defect entries are 
        corrected according to the fixed provided values. 
        More details can be found in https://doi.org/10.1103/PhysRevB.106.134101 .

        The results of the calculation are stored in the `thermodata` property.

        For the chemical potentials, you must provide either:
            -   reservoirs: Dictionary with oxygen partial pressures as keys and dictionary with chemical potential
                as values ({pO2:{'element':chempot}}), or PressureReservoirs object.
        or
            -   precursors + oxygen_ref: Dictionary with {formula:energy} for synthesis precursors and oxygen reference chempot at 0 K.

        Parameters
        ----------
        bulk_dos : (dict or Dos)
            Density of states to integrate. Can be provided as density of states D(E)
            or using effective masses.
            Format for effective masses:
                dict with following keys:
                    - "m_eff_h" : holes effective mass in units of m_e (electron mass)
                    - "m_eff_e" : electrons effective mass in units of m_h          
                    - `band_gap` needs to be provided in args
            Formats for explicit DOS:
                dictionary with following keys:
                    - 'energies' : list or np.array with energy values
                    - 'densitites' : list or np.array with total density values
                    - 'structure' : pymatgen Structure of the material, 
                                    needed for DOS volume and charge normalization.
                or a pymatgen Dos object (Dos, CompleteDos or FermiDos).
        temperature: (float)
            Temperature in K.
        quench_temperature : (float)
            Value of quenching temperature in K.
        quenched_species : (list), optional
            List of defect species to quench. If None all defect species are quenched.The default is None.
        quench_elements : (bool)
            If True the total concentrations of elements at high temperature go in the charge neutrality at low temperature.
            If False the quenched concentrations are the ones of single defect species (e.g. elements are not allowed
            to equilibrate on different sites). The default is False.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names, values are the concentrations. (ex {'Vac_Na':1e20}). 
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
            Must either be a list of dictionaries with {'charge': float, 'conc': float} or a list of SingleDefConc objects. 
        reservoirs: (dict or PressureReservoirs)
            Dictionary with pO2 (float) as keys and chemical potential dictionaly as values ({pO2:{'element':mu_element}})
        precursors: (dict)
            Dictionary with formulas (str) as keys and total energies as values. Chemical potentials are found from the energies of the 
            precursors and the oxygen chempot value (uses the np.linalg.lstsq function). 
            If the system is underdetermined the minimum-norm solution is found.
        oxygen_ref : (float)
            Absolute chempot of oxygen at 0K.
        pressure_range : (tuple)
            Range in which to evaluate the partial pressure.
        npoints : (int)
            Number of data points to interpolate the partial pressure with.
        xtol : (float)
            Tolerance for bisect (scipy) to solve charge neutrality.
        eform_kwargs : (dict)
            Kwargs to pass to `entry.formation_energy`.
        dconc_kwargs : (dict)
            Kwargs to pass to `entry.defect_concentration`.
        kwargs: (dict)
            Kwargs to pass to `plot_pO2_vs_concentrations`.
        """
        from .thermodynamics import DefectThermodynamics
        
        if not reservoirs:
            if not precursors:
                if self.elements != ['O']:
                    raise ValueError('You need to either directly provide reservoirs, or precursors + oxygen chempot reference')
            if not oxygen_ref:
                raise ValueError('You need to provide the oxygen chempot reference when using precursors')

            if self.elements == ['O']:
                reservoirs = get_oxygen_pressure_reservoirs(
                                                            oxygen_ref=oxygen_ref,
                                                            temperature=temperature,
                                                            pressure_range=pressure_range,
                                                            npoints=npoints)
            else:
                reservoirs = get_pressure_reservoirs_from_precursors(
                                                                    precursors=precursors,
                                                                    oxygen_ref=oxygen_ref,
                                                                    temperature=temperature,
                                                                    pressure_range=pressure_range,
                                                                    npoints=npoints)
                
        defects_analysis =  DefectThermodynamics(
                                                defects_analysis=self,
                                                bulk_dos=bulk_dos,
                                                fixed_concentrations=fixed_concentrations,
                                                external_defects=external_defects,
                                                xtol=xtol,
                                                eform_kwargs=eform_kwargs,
                                                dconc_kwargs=dconc_kwargs)
        
        if quench_temperature:
            thermodata = defects_analysis.get_pO2_quenched_thermodata(
                                                                    reservoirs=reservoirs,
                                                                    initial_temperature=temperature,
                                                                    final_temperature=quench_temperature,
                                                                    quenched_species=quenched_species,
                                                                    quench_elements=quench_elements,
                                                                    name='QuenchedBrowerDiagram')
        else:
            thermodata = defects_analysis.get_pO2_thermodata(
                                                            reservoirs=reservoirs,
                                                            temperature=temperature,
                                                            name='BrowerDiagram')
            
        self._thermodata = thermodata
        if 'xlim' not in kwargs.keys():
            kwargs['xlim'] = pressure_range
        plt = plot_pO2_vs_concentrations(thermodata=thermodata,**kwargs)
        
        return plt
    
    
    def plot_doping_diagram(self,
                            variable_defect_specie,
                            concentration_range,
                            chemical_potentials,
                            bulk_dos,
                            temperature,
                            quench_temperature=None,
                            quenched_species=None,
                            quench_elements=False,
                            fixed_concentrations=None,
                            external_defects=[],
                            npoints=50,
                            xtol=1e-10,
                            eform_kwargs={},
                            dconc_kwargs={},
                            **kwargs):
        """
        Calculate defect and carrier concentrations as a function of the concentration of a particular 
        defect species (usually a dopant).

        Parameters
        ----------
        variable_defect_specie : (str)
            Name or element of the variable defect species.
        concentration_range : (tuple or list)
            Logaritmic range of the concentration of the variable species in cm^-3 (ex. [1,20]).
        chemical_potentials : (dict)
            Dictionary containing chemical potentials ({element:chempot}).
        bulk_dos: (dict or Dos)
            Density of states to integrate. Can be provided as density of states D(E)
            or using effective masses.
            Format for effective masses:
                dict with following keys:
                    - "m_eff_h" : holes effective mass in units of m_e (electron mass)
                    - "m_eff_e" : electrons effective mass in units of m_h          
                    - `band_gap` needs to be provided in args
            Formats for explicit DOS:
                dictionary with following keys:
                    - 'energies' : list or np.array with energy values
                    - 'densitites' : list or np.array with total density values
                    - 'structure' : pymatgen Structure of the material, 
                                    needed for DOS volume and charge normalization.
                or a pymatgen Dos object (Dos, CompleteDos or FermiDos).
        temperature: (float)
            Temperature in K.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations. (ex {'Vac_Na':1e20}) 
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
            Must either be a list of dictionaries with {'charge': float, 'conc': float} or a list of SingleDefConc objects. 
        xtol : (float)
            Tolerance for bisect (scipy) to solve charge neutrality.
        npoints : (int), optional
            Number of points to divide concentration range.
        eform_kwargs : (dict)
            Kwargs to pass to `entry.formation_energy`.
        dconc_kwargs : (dict)
            Kwargs to pass to `entry.defect_concentration`.
        kwargs: (dict)
            Kwargs to pass to `plot_variable_species_vs_concentrations`.
        """
        from .thermodynamics import DefectThermodynamics
        
        defects_analysis = DefectThermodynamics(
                                                defects_analysis=self,
                                                bulk_dos=bulk_dos,
                                                fixed_concentrations=fixed_concentrations,
                                                external_defects=external_defects,
                                                xtol=xtol,
                                                eform_kwargs=eform_kwargs,
                                                dconc_kwargs=dconc_kwargs)
        
        if quench_temperature:
            thermodata = defects_analysis.get_variable_species_quenched_thermodata(
                                                        variable_defect_specie=variable_defect_specie,
                                                        concentration_range=concentration_range,
                                                        chemical_potentials=chemical_potentials,
                                                        initial_temperature=temperature,
                                                        final_temperature=quench_temperature,
                                                        quenched_species=quenched_species,
                                                        quench_elements=quench_elements,
                                                        npoints=npoints,
                                                        name='QuenchedDopingDiagram')  
                
        thermodata = defects_analysis.get_variable_species_thermodata(
                                                        variable_defect_specie=variable_defect_specie,
                                                        concentration_range=concentration_range,
                                                        chemical_potentials=chemical_potentials,
                                                        temperature=temperature,
                                                        npoints=npoints,
                                                        name='DopingDiagram')  
        self._thermodata = thermodata
        if 'xlim' not in kwargs.keys():
            kwargs['xlim'] = concentration_range
        plt = plot_variable_species_vs_concentrations(thermodata, **kwargs)

        return plt
        
         
    def plot_formation_energies(self,
                                chemical_potentials,
                                temperature=0,
                                entries=None,
                                xlim=None,
                                ylim=None,
                                title=None,
                                fermi_level=None,
                                grid=True,
                                figsize=(6,6),
                                fontsize=12,
                                show_legend=True,
                                format_legend=True,
                                **eform_kwargs):
        """
        Produce defect Formation energy vs Fermi energy plot.

        Parameters
        -----------        
        chemical_potentials : (dict)
            Dictionary with chemical potentials of the elements {'element':chempot}
        temperature : (float)
            Temperature in K. If no custom formation energy is provided, this arg has no effect.
        entries : (list) 
            List of entries to calculate. If None all entries are considered.
        xlim : (tuple)
            Tuple (min,max) giving the range of the x (fermi energy) axis.
        ylim : (tuple)
            Tuple (min,max) giving the range for the formation energy axis.
        title : (str)
            Title of the figure.
        fermi_level : (float)
            Plot Fermi energy position with a vertical line.
        grid : (bool)
            Show grid.
        figsize : (float or tuple)
            Can be float or tuple for different x and y sizes multiplier to change the figure size.
        fontsize : (float)
            Multiplier to change font size.
        show_legend  : (Bool)
            Show legend.
        format_legend : (bool)
            Get latex-like legend based on the name of defect entries.
        eform_kwargs : (dict)
            Kwargs to pass to `entry.formation_energy`.
                
        --------
        Returns:
            matplotlib object
        """
        entries = entries or self.entries
        kwargs = {
            'entries':entries,
            'chemical_potentials':chemical_potentials,
            'temperature':temperature,
            'vbm':self.vbm,
            'band_gap':self.band_gap,
            'xlim':xlim,
            'ylim':ylim,
            'title':title,
            'fermi_level':fermi_level,
            'grid':grid,
            'figsize':figsize,
            'fontsize':fontsize,
            'show_legend':show_legend,
            'format_legend':format_legend,
            'get_subplot':False,
            'subplot_settings':None,
            'eform_kwargs':eform_kwargs
            }
        if chemical_potentials:
            values = list(chemical_potentials.values())
            if type(values[0]) == dict:
                ncolumns = 2
                nrows = len(values) // ncolumns + len(values) % ncolumns 
                print(len(values),nrows,ncolumns)
                res = chemical_potentials
                idx = 0
                for key,value in res.items():
                    idx += 1
                    kwargs['title'] = key
                    kwargs['chemical_potentials'] = value
                    kwargs['get_subplot'] = True
                    kwargs['subplot_settings'] = (nrows,ncolumns,idx)
                    plt = plot_formation_energies(**kwargs)
            else:
                plt = plot_formation_energies(**kwargs)
        else:
            plt = plot_formation_energies(**kwargs)
        return plt
     
        
    def plot_binding_energies(self,
                            temperature=0,
                            names=None,
                            xlim=None,
                            ylim=None,
                            figsize=(6,6),
                            fontsize=18,
                            format_legend=True,
                            **eform_kwargs):
        """
        Plot binding energies for complex of defects as a function of the fermi level

        Parameters
        ----------
        temperature : (float)
            Temperature in K. If no custom formation energy is provided, this arg has no effect.
        names : (list)
            List of strings with names of DefectEntry. If None, all defect complexes are plotted.
        xlim : (tuple)
            Tuple (min,max) giving the range of the x (fermi energy) axis
        ylim : (tuple)
            Tuple (min,max) giving the range for the formation energy axis
        figsize : (tuple)
            Figure size.
        fontsize : (float)
            Font size.
        format_legend : (bool)
            Bool for getting latex-like legend based on the name of defect entries.

        Returns
        -------
            matplotlib object
        """        
        return plot_binding_energies(entries=self.entries,
                                    vbm=self.vbm,
                                    band_gap=self.band_gap,
                                    temperature=0,
                                    names=names,
                                    xlim=xlim,
                                    ylim=ylim,
                                    figsize=figsize,
                                    fontsize=fontsize,
                                    format_legend=format_legend,
                                    **eform_kwargs)
    
    
    def plot_ctl(self,
                entries=None,
                temperature=0,
                ylim=None,
                figsize=(10,10),
                fontsize=16,
                fermi_level=None,
                format_legend=True,
                get_integers=True,
                **eform_kwargs):
        """
        Plotter for the charge transition levels.

        Parameters
        ----------
        entries : (list)
            List of entries to calculate. If None all entries are considered.
        temperature : (float)
            Temperature in K. If no custom formation energy is provided, this arg has no effect.
        ylim : (tuple)
            y-axis limits.
        figsize : (tuple)
            Figure size.
        fontsize : (float)
            Font size.
        fermi_level : (float)
            Plot Fermi energy position.
        format_legend : (bool)
            Bool for getting latex-like legend based on the name of defect entries.
        get_integers : (bool)
            Get charge transition levels as integers.

        Returns
        -------
            matplotlib object    
        """        
        entries = entries or self.entries
        return plot_charge_transition_levels(entries=entries,
                                            vbm=self.vbm,
                                            band_gap=self.band_gap,
                                            ylim=ylim,
                                            figsize=figsize,
                                            fontsize=fontsize,
                                            fermi_level=fermi_level,
                                            format_legend=format_legend,
                                            get_integers=get_integers,
                                            **eform_kwargs)


    def select_entries(self,
                       entries=None,
                       mode='and',
                       exclude=False,
                       types=None,
                       elements=None,
                       names=None,
                       function=None,
                       **kwargs):
        """
        Find entries based on different criteria. Returns a list of DefectEntry objects.

        Parameters
        ----------
        entries : (list)
            List of defect entries.
        mode : (str), optional
            Filtering mode, possibilities are: 'and' and 'or'. The default is 'and'. 
        exclude : (bool), optional
            Exclude the entries satisfying the criteria instead of selecting them. The default is False.
        types : (list), optional
            Class name of the defect in the entry.
        elements : (list), optional
            List of symbols of elements that need to belong to the defect specie.
            If None this criterion is ignored. The default is None.
        names : (list)
            List of entry names.
        function : (function), optional
            Specific funtion for more complex criteria. The function must take a DefectEntry
            object as argument and return a bool.
        **kwargs : (dict)
            Properties that the jobs need to satisfy. Keys are referred to attributes/methods 
            of the defect entry. To address more than one condition relative to
            the same attribute, use lists or tuples (e.g. charge=[0,1]).

        Returns
        -------
        List of DefectEntry objects.
        """        
        input_entries = entries if entries else self.entries 
        functions = []
        
        if types:
            def ftypes(entry):
                return entry.defect_type in types
            functions.append(ftypes)
        
        if elements:
            def felements(entry):
                for defect in entry.defect:
                    if defect.specie in elements:
                        return True
            functions.append(felements)
                
        if names:
            def fnames(entry):
                return entry.name in names
            functions.append(fnames)
        
        if function:
            functions.append(function)
            
        return select_objects(objects=input_entries,mode=mode,exclude=exclude,
                              functions=functions,**kwargs)
    
    
    def set_formation_energy_functions(self,function,**kwargs):
        """
        Customize the function for formation energy calculations for specific
        defect entries. Pass kwargs to `select_entries` function to choose 
        which entries to apply the function to.
        """
        entries = self.select_entries(**kwargs)
        for entry in entries:
            entry.set_formation_energy_function(function)
        
    def reset_formation_energy_functions(self,reset_all=True,**kwargs):
        """
        Reset formation energy function to default for selected entries, 
        Pass kwargs to `select_entries` function to choose which entries 
        to reset. If reset_all is True all entries are reset to default.
        """
        if reset_all:
            kwargs = {'exclude':True}
        entries = self.select_entries(**kwargs)
        for entry in entries:
            entry.reset_formation_energy_function()

    def set_defect_concentration_functions(self,function,**kwargs):
        """
        Customize the function for defect concentration calculations for specific
        defect entries. Pass kwargs to `select_entries` function to choose 
        which entries to apply the function to.
        """
        entries = self.select_entries(**kwargs)
        for entry in entries:
            entry.set_defect_concentration_function(function)
        
    def reset_defect_concentration_functions(self,reset_all=True,**kwargs):
        """
        Reset defect concentration function to default for selected entries, 
        Pass kwargs to `select_entries` function to choose which entries 
        to reset. If reset_all is True all entries are reset to default.
        """
        if reset_all:
            kwargs = {'exclude':True}
        entries = self.select_entries(**kwargs)
        for entry in entries:
            entry.reset_defect_concentration_function()

    def reset_all_custom_functions(self,reset_all=True,**kwargs):
        """
        Reset all customized functions to default for selected entries, 
        Pass kwargs to `select_entries` function to choose which entries 
        to reset. If reset_all is True all entries are reset to default.
        """
        self.reset_formation_energy_functions(reset_all=reset_all,**kwargs)
        self.reset_defect_concentration_functions(reset_all=reset_all,**kwargs)



    def solve_fermi_level(self,
                        chemical_potentials,
                        bulk_dos,
                        temperature=300,
                        fixed_concentrations=None,
                        external_defects=[],
                        xtol=1e-05,
                        eform_kwargs={},
                        dconc_kwargs={}):
        """
        Solve charge neutrality and get the value of Fermi level at thermodynamic equilibrium.
        
        Parameters
        ----------
        chemical_potentials : (Dict)
            Dictionary of chemical potentials in the format {'element':chempot}.
        bulk_dos : (dict or Dos)
            Density of states to integrate. Can be provided as density of states D(E)
            or using effective masses.
            Format for effective masses:
                dict with following keys:
                    - "m_eff_h" : holes effective mass in units of m_e (electron mass)
                    - "m_eff_e" : electrons effective mass in units of m_h          
                    - `band_gap` needs to be provided in args
            Formats for explicit DOS:
                dictionary with following keys:
                    - 'energies' : list or np.array with energy values
                    - 'densitites' : list or np.array with total density values
                    - 'structure' : pymatgen Structure of the material, 
                                    needed for DOS volume and charge normalization.
                or a pymatgen Dos object (Dos, CompleteDos or FermiDos).
        temperature : (float)
            Temperature in Kelvin.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations (ex {'Vac_A':1e20}) . For more info, 
            read the documentation of the defect_concentrations method.
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
            Must either be a list of dictionaries with {'name': str, 'charge': float, 'conc': float} 
            or a list of SingleDefConc objects. 
        xtol : (float)
            Tolerance for bisect (scipy) to solve charge neutrality.
        eform_kwargs : (dict)
            Kwargs to pass to `entry.formation_energy`.
        dconc_kwargs : (dict)
            Kwargs to pass to `entry.defect_concentration`.

        Returns
        -------
        (float)
            Fermi level satisfying charge neutrality.
        """
                
        def _get_total_q(ef):
            qd_tot = self._get_total_charge(fermi_level=ef,
                                            chemical_potentials=chemical_potentials,
                                            bulk_dos=bulk_dos,
                                            temperature=temperature,
                                            fixed_concentrations=fixed_concentrations,
                                            external_defects=external_defects,
                                            eform_kwargs=eform_kwargs,
                                            dconc_kwargs=dconc_kwargs)
            return qd_tot
        
        root = bisect(_get_total_q, -1, self.band_gap + 1.,xtol=xtol) # set full_output=True for bisect info 
    
        qd_tot = _get_total_q(root)
        if abs(qd_tot) > 1e10:
            warnings.warn(
                    f"Fermi level solver with xtol={xtol} yields high residual charge: "
                    f"({qd_tot: .2e} cm^-3). Check total_charge vs fermi_level behaviour")      

        return root
    
    
    def sort_entries(self,inplace=False,entries=None,features=['name','charge'],reverse=False):
        """
        Sort defect entries with different criteria.

        Parameters
        ----------
        inplace : (bool), optional
            Reset the self.entries attibute with sorted entries. Only works if entries input
            is not given. The default is False.
        entries : (list), optional
            List of defect entries to sort. If None self.entries is used. The default is None.
        features : (list), optional
            List of strings with attribute/method names. The default is ['name','charge'].
        reverse : (bool)
            Reverse order.

        Returns
        -------
        (list)
            List of defect entries.
        """
        inplace = False if entries else inplace
        entries = entries if entries else self.entries
        sorted_entries = sort_objects(objects=entries, features=features)
        if inplace:
            self.entries = sorted_entries
        else:
            return sorted_entries
            
    
    def stable_charges(self,
                       chemical_potentials=None,
                       fermi_level=0,
                       temperature=0,
                       entries=None,
                       **eform_kwargs):
        """
        Creating a dictionary with names of single defect entry as keys and
        as value a tuple (charge,formation_energy) that gives the most stable 
        charge state at the inserted fermi level.

        For each entry, is possible to set a custom function to compute the formation 
        energy using `set_formation_energy_function`. If that is set,
        the new function is used instead of the default. 
        The function args must be the following:
            vbm : (float)
                Valence band maximum of bulk calculation in eV
            chemical_potentials : (dict)
                Chemical potentials of the elements involved in the defect.
            fermi_level : (float)
                Fermi level in eV (with respect to the VBM).
            temperature : (float)
                Temperature in Kelvin.
            kwargs : (dict)
                Additional custom kwargs.
        It is possible to set custom functions collectively using the 
        `set_formation_energy_functions` method in this class.
        To reset functions to default, use `reset_formation_energy_functions`.

        
        Parameters
        ----------
        chemical_potentials: (dict)
            Dictionary of chemical potentials ({element: chempot})   
        fermi_level: (float) 
            Fermi level in eV relative to valence band maximum.
        temperature : (float)
            Temperature in K. If no custom formation energy is provided, this arg has no effect.
        entries : (list)
            List of defect entries to calculate. If None all entries are considered.
        eform_kwargs : (dict)
            Additional custom kwargs to pass to `entry.formation_energy`

        Returns:
            dict in the format {name:(stable charge, formation energy)}
       """
        formation_energies = self.formation_energies(
                                        chemical_potentials=chemical_potentials,
                                        fermi_level=fermi_level,
                                        temperature=temperature,
                                        entries=entries,
                                        **eform_kwargs)        
        stable_charges = {}
        for name in formation_energies:
            emin = 1e40
            for d in formation_energies[name]:
                q = d[0]
                energy = d[1]
                # finding most stable charge state
                if energy < emin:
                    emin = energy
                    q_stable = q
            stable_charges[name] = (q_stable,emin)
            
        return stable_charges
        

    def table(
            self,
            entries=None,
            pretty=False,
            include_bulk=False,
            display=[]):
        """
        Get DataFrame to display entries. 

        Parameters
        ----------
        entries : (list), optional
            Entries to display. If None all entries are displayed. The default is None.
        pretty : (bool), optional
            Optimize DataFrame for prettier visualization.
        include_bulk: (bool), optional
            Include bulk composition and space group for each entry in DataFrame.
        display: (list)
            List of strings with defect entry attributes or method results to display.

        Returns
        -------
        df : 
            DataFrame object.
        """
        if not entries:
            entries = self.entries
        d = {}
        index = []
        table = []
        for e in entries:
            d = {}
            if pretty:
                d['defect'] = e.symbol
            else:
                d['name'] = e.name
            if include_bulk:
                if e.bulk_structure:
                    d['bulk composition'] = e.bulk_structure.composition.formula
                    d['bulk space group'] = e.bulk_structure.get_space_group_info()
                d['bulk_volume'] = e.defect.bulk_volume
            if not pretty:
                d['symbol'] = e.symbol    
                d['delta atoms'] = e.delta_atoms
            d['charge'] = e.charge
            d['multiplicity'] = e.multiplicity
            if pretty:
                if e.corrections:
                    for key,value in e.corrections.items():
                        d[f'{key} correction'] = value
            else:
                d['corrections'] = e.corrections
            if display:
                for feature in display:
                    if isinstance(feature,list):
                        key = feature[0]
                        for k in feature[1:]:
                            key = key + '["%s"]'%k
                    else:
                        key = feature
                    d[key] = get_object_feature(e,feature)
            table.append(d)
        df = pd.DataFrame(table)
        return df
    

    def _get_total_charge(self,
                          fermi_level,
                          chemical_potentials,
                          bulk_dos,
                          temperature=300, 
                          fixed_concentrations=None,
                          external_defects=[],
                          eform_kwargs={},
                          dconc_kwargs={}): 
        """
        Calculate the total charge concentration (defects + holes - electrons) needed to solve charge neutrality.
        Check solve_fermi_level docs
        """
        # defect contribution
        qd_tot = sum([
            d.charge * d.conc
            for d in self.defect_concentrations(
                                        chemical_potentials=chemical_potentials,
                                        temperature=temperature,
                                        fermi_level=fermi_level,
                                        fixed_concentrations=fixed_concentrations,
                                        per_unit_volume=True,
                                        eform_kwargs=eform_kwargs,
                                        **dconc_kwargs)
            ])
        
        # defects with fixed behaviour
        for d_ext in external_defects:
            qd_tot += d_ext['charge'] * d_ext['conc']

        # charge carriers (holes and electrons)
        h, n = get_carrier_concentrations(
                                dos=bulk_dos,
                                fermi_level=fermi_level,
                                temperature=temperature,
                                band_gap=self.band_gap)
        qd_tot += h - n
        return qd_tot



class SingleDefConc(MSONable):
    
    def __init__(self,name,charge,conc):
        """
        Object to store defect concentrations data. It is also subscribtable like a dictionary.

        Parameters
        ----------
        name : (str)
            Name of defect specie.
        charge : (int)
            Charge of defect specie.
        conc : (float)
            Concentration value.
        """
        self.name = name
        self.charge = charge
        self.conc = conc
        
    def __repr__(self):
        s = 'charge=%.1f, conc=%.2e, name=%s' %(self.charge,self.conc,self.name)
        return s
    
    def __print__(self):
        return self.__repr__()
    
    def __getitem__(self,key):
        return getattr(self,key)
    
    def keys(self):
        """
        Returns dictionary-like keys of the object's attributes.
        """
        return ('name', 'charge', 'conc')

    def values(self):
        """
        Returns dictionary-like values of the object's attributes.
        """
        return (self.name, self.charge, self.conc)

    def items(self):
        """
        Returns dictionary-like key-value pairs of the object's attributes.
        """
        return zip(self.keys(), self.values())
        
        
class DefectConcentrations:
    
    def __init__(self,concentrations):
        """
        Class to store sets of defect concentrations (output of concentration calculations with DefectsAnalysis).
        List of SingleDefConc objects. Subscriptable like a list.

        Parameters
        ----------
        concentrations : (list)
            List of SingleDefConc objects.
        """
        for c in concentrations:
            if type(c) == dict:
                c = SingleDefConc(c)
        self.concentrations = concentrations
        self._compute_totals()

    def _compute_totals(self):
        d = {}
        # calculate total concentrations
        for c in self.concentrations:
            gname = c.name # previously was the name without label
            if gname not in d.keys():
                d[gname] = 0
            d[gname] += c.conc    
        self._total = d
      #  store stable concentrations
        conc_stable = []
        for n in self.names:
            concs = self.select_concentrations(name=n)
            cmax = SingleDefConc(name='',conc=None,charge=0) # dummy object
            for c in concs:
                if cmax.conc is None or c.conc >= cmax.conc:
                    cmax = c
            conc_stable.append(concs[concs.index(cmax)])
        self._stable = conc_stable
        return 

    def append(self, item):
        if type(item) == dict:
            item = SingleDefConc.from_dict(item)
        self.concentrations.append(item)
        self._compute_totals()

    def __len__(self):
        return len(self.concentrations)
        
    def __iter__(self):
        return self.concentrations.__iter__()
    
    def __getitem__(self,i):
        return self.concentrations[i]
    
    def __repr__(self):
        return self.__print__()

    def __print__(self):
        return '['+ '\n'.join([c.__print__() for c in self.concentrations])+']'  

    def copy(self):
        return DefectConcentrations(self.concentrations.copy())
   
    def as_dict(self):
        d = [c.as_dict() for c in self.concentrations]
        return d

    @classmethod 
    def from_dict(cls,d):
        dc = [SingleDefConc.from_dict(c) for c in d]
        return cls(dc)
 
    @property
    def elemental(self):
        """
        Dictionary with element (or element vacancy) as keys and total element concentration as values.
        """
        d = {}
        for c in self:
            df = Defect.from_string(c.name)
            for defect in df:
                if defect.specie not in d.keys() or defect.name not in d.keys():
                    if defect.type == 'Vacancy':
                        ekey = defect.name
                        d[ekey] = self.get_element_total(defect.specie,vacancy=True)
                    else:
                        ekey = defect.specie
                        d[ekey] = self.get_element_total(defect.specie,vacancy=False)
        return d
                    
    
    @property
    def elements(self):
        """
        List of element (or vacancies) symbols.
        """
        return self.elemental.keys()

                
    @property
    def names(self):
        """
        List of names of all defect species
        """
        return list(self.total.keys())

    @property
    def stable(self):
        """
        Get concentrations of only stable charge states. List of SingleDefConc objects.
        """
        return DefectConcentrations(self._stable)    

    @property
    def total(self):
        """
        Get total concentrations for every defect specie. Format is a dict ({'name':conc})
        """
        return self._total

    
    def filter_concentrations(self,inplace=False,mode='and',exclude=False,names=None,
                              charges=None,indexes=None,function=None,**kwargs):
        """
        Filter concentrations based on different criteria.

        Parameters
        ----------
        inplace : (bool), optional
            If True update the current object, otherwise returns a new DefectConcentrations object.
            The default is False.
        mode : (str), optional
            Selection mode. available are 'and','or'. The default is 'and'.
        exclude : (bool), optional
            Exclude the entries matching the criteria. The default is False.
        names : (list), optional
            Names of the defect.
        charges : (list), optional
            Charges of the defect.
        indexes : (list), optional
            Indexes of defect in the concentrations list.
        function : (function), optional
            Specific funtion for more complex criteria. The function must take a 
            SingleDefConc object as argument and return a bool.
        **kwargs : (dict)
            Criteria for selection. They need to be attributes of SingleDefConc.

        Returns
        -------
        output_concs : (list)
            Filtered DefectConcentrations object.
        """
        output_concs = self.select_concentrations(mode=mode,exclude=exclude,names=names,
                                                  charges=charges,indexes=indexes,
                                                  function=function,**kwargs)
                    
        if inplace:
            self.concentrations = output_concs
            return
        else:
            return DefectConcentrations(output_concs)
            
        
    def get_element_total(self,element,vacancy=False):
        """
        Get total concentration of every element (summed also across different species)

        Parameters
        ----------
        element : (str)
            Target element.
        vacancy : (bool), optional
            Weather the target element is a vacancy. The default is False.

        Returns
        -------
        eltot : (float)
            Total concentration of target element.
        """
        eltot = 0
        for c in self:
            df = Defect.from_string(c.name)
            for defect in df:
                if element == defect.specie:
                    if vacancy and defect.type=='Vacancy':
                        eltot += c.conc
                    elif vacancy == False and defect.type!='Vacancy':
                        eltot += c.conc
        return eltot
        

    def select_concentrations(self,concentrations=None,mode='and',exclude=False,names=None,
                    charges=None,indexes=None,function=None,**kwargs):
        """
        Select concentrations based on different criteria.

        Parameters
        ----------
        concentrations : (list), optional
            List of SingleDefConc to select from. If None self.concentrations is used.
            The default is None.
        mode : (str), optional
            Selection mode. available are 'and','or'. The default is 'and'.
        exclude : (bool), optional
            Exclude the entries matching the criteria. The default is False.
        names : (list), optional
            Names of the defect.
        charges : (list), optional
            Charges of the defect.
        indexes : (list), optional
            Indexes of defect in the concentrations list.
        function : (function), optional
            Specific funtion for more complex criteria. The function must take a 
            SingleDefConc object as argument and return a bool.
        **kwargs : (dict)
            Criteria for selection. They need to be attributes of SingleDefConc.

        Returns
        -------
        output_concs : (list)
            List of SingleDefConc objects.
        """
        input_concs = concentrations if concentrations else self.concentrations.copy()
        functions = []
        
        if names is not None:
            def fnames(c):
                return c.name in names
            functions.append(fnames)
        
        if charges is not None:
            def fcharges(c):
                return c.charge in charges
            functions.append(fcharges)
        
        if indexes is not None:
            def findexes(c):
                return input_concs.index(c) in indexes
            functions.append(findexes)
        
        if function is not None:
            functions.append(function)
            
        return select_objects(objects=input_concs,mode=mode,exclude=exclude,
                              functions=functions,**kwargs)
    
                    


        
        
        
