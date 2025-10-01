
import numpy as np
from scipy.optimize import bisect
from monty.json import MontyDecoder, MSONable
import pandas as pd
import os
import os.path as op
import json
import copy

from pymatgen.electronic_structure.dos import FermiDos

from .pmg.pmg_dos import FermiDosCarriersInfo
from .chempots.oxygen import get_pressure_reservoirs_from_precursors
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


class DefectsAnalysis:
    
    def __init__(self, entries, vbm, band_gap, sort_entries=True):
        """ 
        Class to compute defect properties starting from single calculations of defects.
        Args:
            entries: (list) 
                A list of DefectEntry objects.
            vbm: (float) 
                Valence Band energy to use for all defect entries.
                NOTE if using band shifting-type correction then this VBM
                should still be that of the GGA calculation
                (the bandedgeshifting_correction accounts for shift's
                contribution to formation energy).
            band_gap : (float)
                Band gap to use for all defect entries.
                NOTE if using band shifting-type correction then this gap
                should still be that of the Hybrid calculation you are shifting to.
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
        Result of thermodynamics calculation (plot_brouwer_diagram or plot_doping_diagram functions).
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
        "@module": self.__class__.__module__,
        "@class": self.__class__.__name__,
        "entries": [entry.as_dict() for entry in self.entries],
        "vbm":self.vbm,
        "band_gap":self.band_gap
            }
        return d


    @classmethod
    def from_dict(cls,d):
        """
        Reconstruct a DefectsAnalysis object from a dict representation created using
        as_dict().

        Args:
            d (dict): dict representation of DefectsAnalysis.

        Returns:
            DefectsAnalysis object
        """
        entries = [MontyDecoder().process_decoded(entry_dict) for entry_dict in d['entries']]        
        vbm = d['vbm']
        band_gap = d['band_gap']
        return cls(entries,vbm,band_gap)


    @staticmethod
    def from_dataframe(df,vbm, band_gap):
        """
        Create DefectsAnalysis object from pandas DataFrame. If df has been 
        exported with include_structure=True, DefectEntry objects will be 
        imported including defect and bulk structures.

        Parameters
        ----------
        df : (DataFrame)
            pandas DataFrame. Needs to be in the following format:
            - rows are defect entries
            - non-optional columns are 'name','charge','multiplicity','energy_diff','bulk_volume'.
                - name : defect name in the format given by Defect objects. 
                - charge : defect charge
                - multiplicity : site multiplicity of the defect
                - energy_diff : difference in energy btw defect and bulk structures
            - If the df was not exported from this class, you need to provide also:
                - bulk_volume : Volume of the bulk cell (unit cell or supercell 
                  relative to your multiplicity in A°3)
        vbm : (float)
            Valence band maximum of bulk structure in eV
        band_gap : (float)
            Band gap of bulk structure in eV
        
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
    def from_file(filename,vbm,band_gap,format=None): 
        """
        Create DefectsAnalysis object from file. Available formats are
        'json','pkl','csv' and 'xlsx'. Check docs in from_dataframe function
        for file formatting requirements.
        """
        if '.json' in filename or format == 'json':
            return DefectsAnalysis.from_json(filename)
        elif '.pkl' in filename or format == 'pkl':
            df = pd.read_pickle(filename)
        elif '.csv' in filename or format == 'csv':
            df = pd.read_csv(filename)
        elif '.xlsx' in filename or format == 'excel':
            df = pd.read_excel(filename)
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
                            get_corrections=False,
                            get_multiplicity=False,
                            get_data=True,
                            band_gap=None,
                            vbm=None,
                            function=None,
                            initial_structure=False,
                            computed_entry_kwargs={},
                            finder_kwargs={}):
        """
        Generate DefectsAnalysis object from VASP directories read with Pymatgen.

        Parameters
        ----------
        path_defects : (str)
            Path of VASP defects calculation.
        path_bulk : (str)
            Path of VASP bulk calculation. 
        get_corrections : (bool)

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
        kwargs : (dict)
            Kwargs to pass to Vasprun.get_computed_entry.

        Returns
        -------
        DefectsAnalysis object.
        """ 
        from pymatgen.io.vasp.outputs import Vasprun
        if band_gap and vbm:
            parse_eigen = False
        else:
            parse_eigen = True
        vasprun_bulk = Vasprun(op.join(path_bulk,'vasprun.xml'),parse_dos=False,parse_eigen=parse_eigen)
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
        
        return DefectsAnalysis(entries=entries,vbm=vbm,band_gap=band_gap)



    def to_dataframe(self,
            entries=None,
            include_structures=False,
            include_data=True,
            properties=[],
            functions={}):
        """
        Export DefectsAnalysis object as DataFrame. 

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
    

    def to_file(self,filename,format='pkl',export_kwargs={},**kwargs):
        """
        Export DefectsAnalysis object to file. 

        Parameters
        ----------
        filename : (str)
            Path of file.
        format : (str)
            Formats available:
            - "pkl" : pickle file. Allows to store structures (default).
            - "json" : DefectsAnalysis object as json file.
            - "csv" : Does not allow to store structures.
            - "excel: Does not allow to store structures.   
        export_kwargs : (dict)
            Kwargs to pass to file exporting function
        kwargs : (dict)
            Kwargs to pass to to_dataframe function. If not provided, 
            'include_strctures' is to True if chosen format is 'pkl', 
            set to False if chosen format is 'csv' or 'excel'.
        """
        if format == 'json':
            self.to_json(path=filename) 
        elif format == 'pkl':
            if 'include_structures' not in kwargs.keys():
                kwargs['include_structures'] = True
            df = self.to_dataframe(**kwargs)
            df.to_pickle(filename,**export_kwargs)
        elif format in ['csv','excel']:
            if 'include_structures' not in kwargs.keys():
                kwargs['include_structures'] = False
            df = self.to_dataframe(**kwargs)
            if format == 'csv':
                df.to_csv(filename,**export_kwargs)
            elif format == 'excel':
                df.to_excel(filename,**export_kwargs)

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
        d = self.as_dict(**kwargs)
        if path == '':
            path = op.join(self.path,self.name+'.json')
        if path:
            with open(path,'w') as file:
                json.dump(d,file)
            return
        else:
            return d.__str__() 
        



    def binding_energy(self,name,fermi_level=0):
        """
        Args:
            name (string): name of defect complex as assigned in defect entry object
            fermi_level (float): Fermi level at which binding energy needs to be computed
        Returns:
            binding_energy (float)
        """
        stable_charges = self.stable_charges(None,fermi_level)
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
        given a fixed Fermi level
        Args:
            bulk_dos: bulk system dos (pymatgen Dos object)
            temperature: Temperature to equilibrate fermi energies for
            fermi_level: (float) is fermi level relative to valence band maximum
                Default efermi = 0 = VBM energy         
        Returns:
            h,n in absolute values
        """
        return get_carrier_concentrations(dos=bulk_dos,fermi_level=fermi_level,temperature=temperature)


    def charge_transition_levels(self, energy_range=None,entries=None,get_integers=True):
        """
        Computes charge transition levels for all defect entries
        Args:
            energy_range (Tuple): Energy range to evaluate the charge transition levels in
                                  default to (-0.5, Eg + 0.5)    
            entries:
                List of entries to calculate. If None all entries are considered.
            get_integers : (bool)
                Save charges as integers. More convenient for plotting.
        Returns:
            Dictionary with defect name and list of tuples for charge transition levels:
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
        previous_stable_charges = self.stable_charges(None,fermi_level=energy_range[0],entries=entries)
        charge_transition_levels = {name:[] for name in previous_stable_charges}
        for i in range(0,len(e)):
            stable_charges = self.stable_charges(None, fermi_level = e[i],entries=entries)
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
    

    def defect_concentrations(self, chemical_potentials, temperature=300, fermi_level=0.,
                              fixed_concentrations=None,per_unit_volume=True):
        """
        Give list of all concentrations at specified efermi.
        If fixed_concentrations is provided the concentration of defect entries are 
        corrected according to the fixed provided values. More details can be found in 
        https://doi.org/10.1103/PhysRevB.106.134101 .
            
        Parameters
        ----------
        chemical_potentials: (dict)
            Dictionary of chemical potentials ({element: chempot})   
        temperature: (float) 
            Temperature to produce concentrations at.
        fermi_level: (float) 
            Fermi level relative to valence band maximum. The default is 0. 
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys can be simple element strings
            (or vacancies of elements in the format 'Vac_{el}') if that element needs to be 
            fixed across all defect species, alternatively defect entry names can be used as well 
            to target specific defect entries. The values are the concentrations.
        per_unit_volume: (bool)
            Get concentrations in cm^-3. If False they are per unit cell. The default is True.
        
        Returns:
        --------
            list of dictionaries of defect concentrations
        """
        concentrations = []
        if fixed_concentrations:
            dc = self.defect_concentrations(chemical_potentials,temperature,fermi_level,
                                            fixed_concentrations=None,per_unit_volume=per_unit_volume)
            frozen = fixed_concentrations 

        for e in self.entries:
            nsites = e.defect.site_concentration_in_cm3 if per_unit_volume else e.multiplicity
            # frozen defects approach
            if fixed_concentrations:
                c = e.defect_concentration(self.vbm, chemical_potentials,temperature,fermi_level,per_unit_volume)
                corr = self._get_frozen_correction(e,frozen,dc)
                c = c * corr
                defconc = SingleDefConc(name=e.name,charge=e.charge,conc=c)
                concentrations.append(defconc)     
            
            else:
                c = e.defect_concentration(self.vbm, chemical_potentials,temperature,fermi_level,per_unit_volume)
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

    
    def filter_entries(self,inplace=False,entries=None,mode='and',exclude=False,types=None,
                       elements=None,names=None,function=None,**kwargs):
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


    def formation_energies(self,chemical_potentials,fermi_level=0,entries=None):
        """
        Creating a dictionary with names of single defect entry as keys and
        a list of tuples (charge,formation_energy) as values
        Every defect name identifies a type of defect
        Args:
            chemical_potentials:
                a dictionnary of {Element:value} giving the chemical
                potential of each element
            fermi_level:
                Value of fermi level to use for calculating formation energies
            entries:
                List of entries to calculate. If None all entries are considered.
        Returns:
            {name: [(charge,formation energy)] }
        """
        formation_energies = {}
        entries = entries if entries else self.entries
        for entry in entries:
            name = entry.name
            charge = entry.charge
            eform = entry.formation_energy(self.vbm,chemical_potentials,fermi_level=fermi_level)
            if name in formation_energies:
                formation_energies[name].append((charge,eform))
            else:
                formation_energies[name] = []
                formation_energies[name].append((charge,eform))
        
        return formation_energies
        
           
    def get_charge_transition_level(self,name,q1,q2):
        """
        Args:
            chemical_potentials: 
                a dictionnary of {element:champot} giving the chemical
                potential of each element
            name: (str) Name of the defect system to get the charge transition levels for
            q1, q2 (int) : Charges for whcih you want to evaluate the transition levels
        Returns: Charge transition level (float)
        """
        chemical_potentials = None
        formation_energies = self.formation_energies(chemical_potentials)
        
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
                            fixed_concentrations=None,
                            external_defects=[],
                            reservoirs = None,
                            precursors = None,
                            oxygen_ref = None,
                            pressure_range=(1e-20,1e10),
                            npoints = 50,
                            xtol=1e-05,
                            **kwargs):
        """
        Plot Brouwer diagram (defect concentrations vs oxygen partial pressure). Wrapper function for DefectThermodynamics
        and ThermodynamicsPlotter, if you need more control use the classes individually. 

        For the chemical potentials, you must provide either:
            -   reservoirs: Dictionary with oxygen partial pressures as keys and dictionary with chemical potential
                as values ({pO2:{'element':chempot}}), or PressureReservoirs object.
        or
            -   precursors + oxygen_ref: Dictionary with {formula:energy} for synthesis precursors and oxygen reference chempot at 0 K.

        Parameters
        ----------
        bulk_dos : (dict or Dos)
            Density of states to integrate. 
            Can either be a dictionary with following keys:
                - 'energies' : list or np.array with energy values
                - 'densitites' : list or np.array with total density values
                - 'structure' : pymatgen Structure of the material, needed for DOS volume normalization.
            or a pymatgen Dos object (Dos, CompleteDos or FermiDos).
        temperature: (float)
            Temperature in K.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names, values are the concentrations. (ex {'Vac_Na':1e20}). 
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
            Must either be a list of dictionaries with {'charge': float, 'conc': float} 
            or a list of SingleDefConc objects. 
        reservoirs: (dict or PressureReservoirs)
            Dictionary with pO2 as keys and chemical potential dictionaly as values ({pO2:{'element':mu_element}})
        precursors: (dict)
            Dictionary with formulas (str) as keys and total energies as values. Chemical potentials are found from the energies of the 
            precursors and the oxygen chempot value (uses the np.linalg.lstsq function). 
            If the system is underdetermined the minimum-norm solution is found.
        oxygen_ref : (float)
            Absolute chempot of oxygen at 0K.
        pressure_range : (tuple)
            Range in which to evaluate the partial pressure. The default is from 1e-20 to 1e10.
        npoints : (int)
            Number of data points to interpolate the partial pressure with. The default is 50.
        xtol: Tolerance for bisect (scipy) to solve charge neutrality. The default is 1e-05.
        kwargs: (dict)
            Kwargs to pass to plot_pO2_vs_concentrations function.
        """
        from .thermodynamics import DefectThermodynamics
        
        if not reservoirs:
            if not precursors:
                raise ValueError('You need to either directly provide reservoirs, or precursors + oxygen chempot reference')
            if not oxygen_ref:
                raise ValueError('You need to provide the oxygen chempot reference when using precursors')
            reservoirs = get_pressure_reservoirs_from_precursors(precursors=precursors,
                                                                 oxygen_ref=oxygen_ref,
                                                                 temperature=temperature,
                                                                 pressure_range=pressure_range,
                                                                 npoints=npoints)
            
        defects_analysis =  DefectThermodynamics(defects_analysis=self,
                                                bulk_dos=bulk_dos,
                                                fixed_concentrations=fixed_concentrations,
                                                external_defects=external_defects,
                                                xtol=xtol)
        thermodata = defects_analysis.get_pO2_thermodata(reservoirs=reservoirs,
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
                            fixed_concentrations=None,
                            external_defects=[],
                            xtol=1e-05,
                            npoints=50,
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
        chemical_potentials : (Chempots)
            Chempots object containing chemical potentials.
        bulk_dos: (Dos)
            Density of state of bulk material (Pymatgen Dos object).
        temperature: (float)
            Temperature in K.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations. (ex {'Vac_Na':1e20}) 
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
            Must either be a list of dictionaries with {'charge': float, 'conc': float} 
            or a list of SingleDefConc objects. 
        xtol: Tolerance for bisect (scipy) to solve charge neutrality. The default is 1e-05.
        npoints : (int), optional
            Number of points to divide concentration range. The default is 50.
        """
        from .thermodynamics import DefectThermodynamics
        
        defects_analysis = DefectThermodynamics(defects_analysis=self,
                                          bulk_dos=bulk_dos,
                                          fixed_concentrations=fixed_concentrations,
                                          external_defects=external_defects,
                                          xtol=xtol)
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
                                entries=None,
                                xlim=None,
                                ylim=None,
                                title=None,
                                fermi_level=None,
                                grid=True,
                                figsize=(6,6),
                                fontsize=12,
                                show_legend=True,
                                format_legend=True):
        """
        Produce defect Formation energy vs Fermi energy plot.

        -----------
        Parameters:

            chemical_potentials : (dict)
                Dictionary with chemical potentials of the elements {'element':chempot}
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
                
        --------
        Returns:
            matplotlib object
        """
        entries = entries or self.entries
        kwargs = {
            'entries':entries,
            'chemical_potentials':chemical_potentials,
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
            'subplot_settings':None
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
        return plt
     
        
    def plot_binding_energies(self, names=None,xlim=None,ylim=None,figsize=(6,6),fontsize=18,format_legend=True):
        """
        Plot binding energies for complex of defects as a function of the fermi level
        Args:
            names: 
                List of strings with names of DefectEntry. If None all DefectEntry
                objects are plotted.
            xlim:
                Tuple (min,max) giving the range of the x (fermi energy) axis
            ylim:
                Tuple (min,max) giving the range for the formation energy axis
            size:
                Float multiplier to change plot size
            format_legend:
                Bool for getting latex-like legend based on the name of defect entries
        Returns:
            matplotlib object
        """        
        return plot_binding_energies(entries=self.entries,
                                     vbm=self.vbm,
                                     band_gap=self.band_gap,
                                     names=names,
                                     xlim=xlim,
                                     ylim=ylim,
                                     figsize=figsize,
                                     fontsize=fontsize,
                                     format_legend=format_legend)
    
    
    def plot_ctl(self,
                 entries=None,
                 ylim=None,
                 figsize=(10,10),
                 fontsize=16,
                 fermi_level=None,
                 format_legend=True,
                 get_integers=True):
        """
        Plotter for the charge transition levels
        Args:
            entries: List of entries to calculate. If None all entries are considered.
            ylim (Tuple): y-axis limits.
            fermi_level (float) : float to plot Fermi energy position.
            size (float) : Float multiplier for plot size.
            format_legend (bool): Bool for getting latex-like legend based on the name of defect entries .
            get_integers (bool): Get CTLs as integers.
        Returns:
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
                                             get_integers=get_integers)


    def select_entries(self,entries=None,mode='and',exclude=False,types=None,elements=None,
                       names=None,function=None,**kwargs):
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
    
    
    def solve_fermi_level(self,chemical_potentials,bulk_dos,temperature=300,
                                    fixed_concentrations=None,external_defects=[],xtol=1e-05):
        """
        Solve charge neutrality and get the value of Fermi level at thermodynamic equilibrium.
        
        Parameters
        ----------
        chemical_potentials : (Dict)
            Dictionary of chemical potentials in the format {'element':chempot}.
        bulk_dos : (dict or Dos)
            Density of states to integrate. 
            Can either be a dictionary with following keys:
                - 'energies' : list or np.array with energy values
                - 'densitites' : list or np.array with total density values
                - 'structure' : pymatgen Structure of the material, needed for DOS volume normalization.
            or a pymatgen Dos object (Dos, CompleteDos or FermiDos).
        temperature : (float)
            Temperature to equilibrate the system to. The default is 300.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations (ex {'Vac_A':1e20}) . For more info, 
            read the documentation of the defect_concentrations method.
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
            Must either be a list of dictionaries with {'name': str, 'charge': float, 'conc': float} 
            or a list of SingleDefConc objects. 
        xtol: Tolerance for bisect (scipy) to solve charge neutrality. The default is 1e-05.

        Returns
        -------
        (float)
            Fermi level dictated by charge neutrality.
        """
                
        def _get_total_q(ef):
            qd_tot = self._get_total_charge(fermi_level=ef,
                                            chemical_potentials=chemical_potentials,
                                            bulk_dos=bulk_dos,
                                            temperature=temperature,
                                            fixed_concentrations=fixed_concentrations,
                                            external_defects=external_defects)
            return qd_tot
                       
        return bisect(_get_total_q, -1., self.band_gap + 1.,xtol=xtol)
    
    
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
            
    
    def stable_charges(self,chemical_potentials,fermi_level=0,entries=None):
        """
        Creating a dictionary with names of single defect entry as keys and
        as value a tuple (charge,formation_energy) that gives the most stable 
        charge state at the inserted fermi level.
        Every defect name identifies a type of defect
        Args:
            chemical_potentials:
                a dictionnary of {Element:value} giving the chemical
                potential of each element
            fermi_level:
                Value of fermi level to use for calculating formation energies 
            entries:
                List of entries to calculate. If None all entries are considered.
        Returns:
            {name:(stable charge, formation energy)}
       """
        formation_energies = self.formation_energies(chemical_potentials,fermi_level,entries)        
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
        

    def table(self,entries=None,pretty=False,include_bulk=False,display=[]):
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
            symbol = e.symbol
            if pretty:
                index.append(symbol)
            else:
                index.append(e.name)
            d = {}
            if include_bulk:
                d['bulk composition'] = e.bulk_structure.composition.formula
                d['bulk space group'] = e.bulk_structure.get_space_group_info()
            if not pretty:
                d['symbol'] = symbol    
                d['delta atoms'] = e.delta_atoms
            d['charge'] = e.charge
            d['multiplicity'] = e.multiplicity
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
        df = pd.DataFrame(table,index=index)
        if pretty:
            df.index.name = 'symbol'
        else:
            df.index.name = 'name'
        return df
    

    def _get_total_charge(self,
                          fermi_level,
                          chemical_potentials,
                          bulk_dos,
                          temperature=300, 
                          fixed_concentrations=None,
                          external_defects=[]): 
        """
        Calculate the total charge concentration (defects + holes - electrons) needed to solve charge neutrality.
        Check solve_fermi_level docs
        """
        qd_tot = sum([
            d.charge * d.conc
            for d in self.defect_concentrations(chemical_potentials,temperature,fermi_level,
                                                fixed_concentrations,per_unit_volume=True)])
        for d_ext in external_defects:
            qd_tot += d_ext['charge'] * d_ext['conc']

        h, n = get_carrier_concentrations(dos=bulk_dos,fermi_level=fermi_level,temperature=temperature,band_gap=self.band_gap)
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
        return self._stable    

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
    
                    


        
        
        
