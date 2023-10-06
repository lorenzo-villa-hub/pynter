#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:48:55 2021

@author: villa
"""
import warnings
from monty.json import MSONable
from monty.json import jsanitize
import numpy as np

from pynter.defects.analysis import DefectConcentrations
import copy
import os.path as op
import json
import os


class Conductivity:
    """
    Class that handles conductivity calculations.
    """
    def __init__(self,mobilities):
        """
        Parameters
        ----------
        mobilities : (dict)
            Dictionary with mobility values for the defect species in [cm^2/(V*s)].
            Keys must contain "electrons", "holes" and the defect species names.
        """
        self.mobilities = mobilities
        
        
    def get_conductivity(self,carrier_concentrations,defect_concentrations):
        """
        Calculate conductivity from the concentrations of electrons, holes and defects and their mobilities.
        
        Parameters
        ----------
        carrier_concentrations : (list)
            List of tuples with intrinsic carriers concentrations (holes,electrons) in cm^-3.
        defect_concentrations : (list)
            Defect concentrations in the same format as the output of DefectsAnalysis in cm^-3. 

        Returns
        -------
        sigma : (float)
            Conductivity in S/m.

        """
        e = 1.60217662e-19
        mob = self.mobilities
        cc = carrier_concentrations
        sigma_el = e * (mob['holes']*cc[0]+ mob['electrons']*cc[1])
        dc = defect_concentrations
        sigma_ionic = 0
        for d in dc:
            dname = d['name']
            if dname in mob.keys():
                sigma_ionic += mob[dname] * d['conc'] * abs(d['charge']) * e
        sigma = sigma_el + sigma_ionic * 1e02 # conversion from S/cm to S/m
        
        return sigma


    def get_conductivities_from_thermodata(self,thermodata):
        """
        Calculate conductivities from defect and carrier concentrations in
        thermodynamic data. Only the defect species which are present in the
        mobility dict are considered for the calculation of the conductivity.

        Parameters
        ----------
        thermodata: (ThermoData)
            ThermoData object that contains the thermodynamic data.

        Returns
        -------
        conductivities: (list)
            List with conductivity values in S/m.
        """
        conductivities = []
        for i in range(0,len(thermodata.carrier_concentrations)):
            carrier_concentrations = thermodata.carrier_concentrations[i]
            dc = thermodata.defect_concentrations[i] 
            defect_concentrations = []
            for d in dc:
                if d['name'] in self.mobilities.keys():
                    defect_concentrations.append(d)
            conductivity = self.get_conductivity(
                                                carrier_concentrations=carrier_concentrations,
                                                defect_concentrations=defect_concentrations)
            conductivities.append(conductivity)
        
        return conductivities 



class DefectThermodynamics:
    """
    Class that handles the analysis of defect equilibrium.
    """
    
    def __init__(self,defects_analysis,bulk_dos,fixed_concentrations=None,
                 external_defects=[],xtol=1e-05):
        """
        Parameters
        ----------
        defects_analysis :
            DefectsAnalysis object.
        bulk_dos : 
            Pymatgen Dos object.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations. (ex {'Vac_Na':1e20}) 
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
        xtol: Tolerance for bisect (scipy) to solve charge neutrality. The default is 1e-05.
        """
        self.da = defects_analysis
        self.bulk_dos = bulk_dos
        self.fixed_concentrations = fixed_concentrations if fixed_concentrations else None
        self.external_defects = external_defects if external_defects else []
        self.xtol = xtol


    def get_pO2_thermodata(self,reservoirs,temperature=None,name=None):
        """
        Calculate defect and carrier concentrations as a function of the oxygen partial pressure.

        Parameters
        ----------
        reservoirs : (dict, Reservoirs or PressureReservoirs)
            Object with partial pressure values as keys and chempots dictionary as values.
        temperature : (float), optional
            Temperature in Kelvin. If None reservoirs.temperature is used. The default is None.
        name : (str), optional
            Name to assign to ThermoData.

        Returns
        -------
        thermodata : (ThermoData)
            ThermoData object that contains the thermodynamic data:
                partial_pressures : (list)
                    List of partial pressure values.
                defect_concentrations : (list)
                    List of DefectConcentrations objects
                carrier_concentrations : (list)
                    List of tuples with intrinsic carriers concentrations (holes,electrons).
                fermi_levels : (list)
                    list of Fermi level values
        """
        res = reservoirs
        if temperature:
            T = temperature
        elif hasattr(res,'temperature'):
            T = res.temperature
        else:
            raise ValueError('Temperature needs to be provided or to be present ad attribute in PressureReservoirs object')
        partial_pressures = list(res.keys())
        defect_concentrations = []
        carrier_concentrations = []
        fermi_levels=[]

        for r,mu in res.items():
            single_thermodata = self.get_single_point_thermodata(chemical_potentials=mu, 
                                                                 temperature=T)

            defect_concentrations.append(single_thermodata['defect_concentrations'])
            carrier_concentrations.append(single_thermodata['carrier_concentrations'])
            fermi_levels.append(single_thermodata['fermi_levels'])
        
        thermodata = {}
        thermodata['partial_pressures'] = partial_pressures 
        thermodata['defect_concentrations'] = defect_concentrations 
        thermodata['carrier_concentrations'] = carrier_concentrations
        thermodata['fermi_levels'] = fermi_levels
        
        thermodata = ThermoData(thermodata,temperature=temperature,name=name)
        
        return thermodata


    def get_pO2_quenched_thermodata(self,reservoirs,initial_temperature,final_temperature,
                                  quenched_species=None,quench_elements=False,name=None):
        """
        Calculate defect and carrier concentrations as a function of oxygen partial pressure with quenched defects. 
        It is possible to select which defect species to quench and which ones are free to equilibrate.
        Frozen defects in the inputs of the class are still considered.

        Parameters
        ----------
        reservoirs : (dict, Reservoirs or PressureReservoirs)
            Object with partial pressure values as keys and chempots dictionary as values.
        initial_temperature : (float)
            Value of initial temperature (K).
        final_temperature : (float)
            Value of final temperature (K).
        quenched_species : (list), optional
            List of defect species to quench. If None all defect species are quenched.The default is None.
        quench_elements : (bool)
            If True the total concentrations of elements at high temperature go in the charge neutrality at low temperature.
            If False the quenched concentrations are the ones of single defect species (e.g. elements are not allowed
            to equilibrate on different sites). The default is False.
        name : (str), optional
            Name to assign to ThermoData.
            
        Returns
        -------
        thermodata : (ThermoData)
            ThermoData object that contains the thermodynamic data:
                partial_pressures : (list)
                    List of partial pressure values.
                defect_concentrations : (list)
                    List of DefectConcentrations objects
                carrier_concentrations : (list)
                    List of tuples with intrinsic carriers concentrations (holes,electrons).
                fermi_levels : (list)
                    list of Fermi level values
        """        
        res = reservoirs
        if hasattr(res,'temperature'):
            if res.temperature != initial_temperature:
                warnings.warn('PressureReservoirs temperature is not set to the initial quenching temperature',UserWarning)
        partial_pressures = list(res.keys())
        fermi_levels = []
        defect_concentrations = []
        carrier_concentrations = []
                
        for r,mu in res.items():
            single_quenched_thermodata = self.get_single_point_quenched_thermodata(
                                                        chemical_potentials=mu,
                                                        initial_temperature=initial_temperature,
                                                        final_temperature=final_temperature,
                                                        quenched_species=quenched_species,
                                                        quench_elements=quench_elements)
            
            carrier_concentrations.append(single_quenched_thermodata['carrier_concentrations'])
            defect_concentrations.append(single_quenched_thermodata['defect_concentrations'])
            fermi_levels.append(single_quenched_thermodata['fermi_levels'])
 
        thermodata = {}
        thermodata['partial_pressures'] = partial_pressures
        thermodata['fermi_levels'] = fermi_levels
        thermodata['defect_concentrations'] = defect_concentrations
        thermodata['carrier_concentrations'] = carrier_concentrations
        
        thermodata = ThermoData(thermodata,temperature=(initial_temperature,final_temperature),name=name)
        
        return thermodata


    def get_single_point_thermodata(self,chemical_potentials,temperature,
                                    fixed_concentrations=None,external_defects=None,
                                    name=None):
        """
        Compute carrier concentrations, defect concentrations and Fermi level for 
        a single set of chemical potentials.

        Parameters
        ----------
        chemical_potentials : (Chempots)
            Chempots object containing chemical potentials.
        temperature : (int)
            Temperature.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations. (ex {'Vac_Na':1e20}) 
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
        name : (str)
            Label for ThermoData. The default is None.

        Returns
        -------
        thermodata : (ThermoData)
            ThermoData object that contains the thermodynamic data:
                defect_concentrations : (DefectConcentrations)
                    DefectConcentrations object.
                carrier_concentrations : (tuple)
                    Tuple with intrinsic carriers concentrations (holes,electrons).
                fermi_levels : (float)
                    Fermi level value
        """
        dos = self.bulk_dos
        fixed_df = fixed_concentrations if fixed_concentrations else self.fixed_concentrations
        ext_df = self.external_defects if external_defects else self.external_defects
        
        fermi_level = self.da.solve_fermi_level(chemical_potentials=chemical_potentials,
                                                bulk_dos=dos,temperature=temperature,
                                                fixed_concentrations=fixed_df,
                                                external_defects=ext_df,xtol=self.xtol)
        
        carrier_concentrations = self.da.carrier_concentrations(bulk_dos=dos,temperature=temperature,
                                                  fermi_level=fermi_level)
        
        defect_concentrations = self.da.defect_concentrations(chemical_potentials=chemical_potentials,
                                                              temperature=temperature,fermi_level=fermi_level,
                                                              fixed_concentrations=fixed_df)
        
        thermodata = {'carrier_concentrations':carrier_concentrations,
                      'defect_concentrations':defect_concentrations,
                      'fermi_levels':fermi_level}
        
        return ThermoData(thermodata,temperature=temperature,name=name)
    

    def get_single_point_quenched_thermodata(self,chemical_potentials,
                                    initial_temperature,final_temperature,
                                    quenched_species=None, quench_elements=False,
                                    fixed_concentrations=None,external_defects=None,
                                    name=None):
        """
        Compute carrier concentrations, defect concentrations and Fermi level for 
        a single set of chemical potentials.

        Parameters
        ----------
        chemical_potentials : (Chempots)
            Chempots object containing chemical potentials.
        initial_temperature : (float)
            Value of initial temperature (K).
        final_temperature : (float)
            Value of final temperature (K).
        quenched_species : (list), optional
            List of defect species to quench. If None all defect species are quenched.The default is None.
        quench_elements : (bool)
            If True the total concentrations of elements at high temperature go in the charge neutrality at low temperature.
            If False the quenched concentrations are the ones of single defect species (e.g. elements are not allowed
            to equilibrate on different sites). The default is False.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations. (ex {'Vac_Na':1e20}) 
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
        name : (str)
            Label for ThermoData. The default is None.

        Returns
        -------
        thermodata : (ThermoData)
            ThermoData object that contains the thermodynamic data:
                defect_concentrations : (DefectConcentrations)
                    DefectConcentrations object.
                carrier_concentrations : (tuple)
                    Tuple with intrinsic carriers concentrations (holes,electrons).
                fermi_levels : (float)
                    Fermi level value
        """
        fixed_df = fixed_concentrations if fixed_concentrations else self.fixed_concentrations
        ext_df = external_defects if external_defects else self.external_defects
        
        single_thermodata = self.get_single_point_thermodata(
                                                    chemical_potentials=chemical_potentials,
                                                    temperature=initial_temperature,
                                                    fixed_concentrations=fixed_df,
                                                    external_defects=ext_df
                                                    )

        if quench_elements:
            c1 = single_thermodata['defect_concentrations'].elemental
        else:
            c1 = single_thermodata['defect_concentrations'].total

        if quenched_species is None:
            quenched_concentrations = c1.copy()
        else:
            if fixed_df:
                quenched_concentrations = copy.deepcopy(fixed_df)
            else:
                quenched_concentrations = {}
            for k in quenched_species:
                quenched_concentrations[k] = c1[k]
        
        single_quenched_thermodata = self.get_single_point_thermodata(
                                                    chemical_potentials=chemical_potentials,
                                                    temperature=final_temperature,
                                                    fixed_concentrations=quenched_concentrations,
                                                    external_defects=ext_df
                                                    )        
        return single_quenched_thermodata


    def get_variable_species_thermodata(self,variable_defect_specie,concentration_range,
                                        chemical_potentials,temperature,npoints=50,name=None):
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
        temperature : (float)
            Temperature.
        npoints : (int), optional
            Number of points to divide concentration range. The default is 50.
        name : (str), optional
            Label for ThermoData. The default is None.

        Returns
        -------
        thermodata : (ThermoData)
            ThermoData object that contains the thermodynamic data:
                variable_defect_specie : (str)
                    Name of variable defect species.
                variable_concentrations : (list)
                    List of concentrations of variable species. 
                defect_concentrations : (list or dict)
                    Defect concentrations in the same format as the output of DefectsAnalysis. 
                carrier_concentrations : (list)
                    List of tuples with intrinsic carriers concentrations (holes,electrons).
                fermi_levels : (list)
                    List of Fermi level values.
        """
        carrier_concentrations = []
        defect_concentrations = []
        fermi_levels = []
        
        concentrations = np.logspace(start=concentration_range[0],stop=concentration_range[1],num=npoints)
        fixed_df = self.fixed_concentrations.copy() if self.fixed_concentrations else {}
        for c in concentrations:
            fixed_df.update({variable_defect_specie:c})
            
            single_thermodata = self.get_single_point_thermodata(
                                            chemical_potentials=chemical_potentials,
                                            temperature=temperature,fixed_concentrations=fixed_df
                                            )
            
            defect_concentrations.append(single_thermodata['defect_concentrations'])
            carrier_concentrations.append(single_thermodata['carrier_concentrations'])
            fermi_levels.append(single_thermodata['fermi_levels'])
            
        data = {}
        data['variable_defect_specie'] = variable_defect_specie
        data['variable_concentrations'] = concentrations
        data['defect_concentrations'] = defect_concentrations
        data['carrier_concentrations'] = carrier_concentrations
        data['fermi_levels'] = fermi_levels
        
        thermodata = ThermoData(data,temperature=temperature,name=name)
        
        return thermodata


    def get_variable_species_quenched_thermodata(self,
                                  variable_defect_specie,concentration_range,
                                  chemical_potentials,initial_temperature,final_temperature,
                                  quenched_species=None,quench_elements=False,npoints=50,name=None):
        """
        Calculate defect and carrier concentrations as a function of the concentration of a particular 
        defect species (usually a dopant) with quenched defects.
        It is possible to select which defect species to quench and which ones are free to equilibrate.
        Frozen defects in the inputs of the class are still considered.

        Parameters
        ----------
        variable_defect_specie : (str)
            Name or element of the variable defect species.
        concentration_range : (tuple or list)
            Logaritmic range of the concentration of the variable species in cm^-3 (ex. [1,20]).
        chemical_potentials : (Chempots)
            Chempots object containing chemical potentials.
        initial_temperature : (float)
            Value of initial temperature (K).
        final_temperature : (float)
            Value of final temperature (K).
        quenched_species : (list), optional
            List of defect species to quench. If None all defect species are quenched.The default is None.
        quench_elements : (bool)
            If True the total concentrations of elements at high temperature go in the charge neutrality at low temperature.
            If False the quenched concentrations are the ones of single defect species (e.g. elements are not allowed
            to equilibrate on different sites). The default is False.
        npoints : (int), optional
            Number of points to divide concentration range. The default is 50.
        name : (str), optional
            Name to assign to ThermoData.

        Returns
        -------
        thermodata : (ThermoData)
            ThermoData object that contains the thermodynamic data:
                variable_defect_specie : (str)
                    Name of variable defect species.
                variable_concentrations : (list)
                    List of concentrations of variable species. 
                defect_concentrations : (list or dict)
                    Defect concentrations in the same format as the output of DefectsAnalysis. 
                carrier_concentrations : (list)
                    List of tuples with intrinsic carriers concentrations (holes,electrons).
                fermi_levels : (list)
                    List of Fermi level values.
        """
        concentrations = np.logspace(start=concentration_range[0],stop=concentration_range[1],num=npoints)
        fermi_levels = []
        defect_concentrations = []
        carrier_concentrations = []
            
        fixed_df = self.fixed_concentrations.copy() if self.fixed_concentrations else {}
        for c in concentrations:
            fixed_df.update({variable_defect_specie:c})
            
            single_quenched_thermodata = self.get_single_point_quenched_thermodata(
                                                    chemical_potentials=chemical_potentials,
                                                    initial_temperature=initial_temperature,
                                                    final_temperature=final_temperature,
                                                    quenched_species=quenched_species,
                                                    quench_elements=quench_elements,
                                                    fixed_concentrations=fixed_df
                                                    )
            
            defect_concentrations.append(single_quenched_thermodata['defect_concentrations'])
            carrier_concentrations.append(single_quenched_thermodata['carrier_concentrations'])
            fermi_levels.append(single_quenched_thermodata['fermi_levels'])
                
            
        data = {}
        data['variable_defect_specie'] = variable_defect_specie
        data['variable_concentrations'] = concentrations
        data['fermi_levels'] = fermi_levels
        data['defect_concentrations'] = defect_concentrations
        data['carrier_concentrations'] = carrier_concentrations
        
        thermodata = ThermoData(data,temperature=(initial_temperature,final_temperature),name=name)
        
        return thermodata



class ThermoData(MSONable):
    
    "Class to handle defect thermodynamics data"
    
    def __init__(self,thermodata,temperature=None,name=None):
        """
        Class that handles dict of defect thermodynamics data.
        
        Parameters
        ----------
        thermodata : (dict)
            Dict that contains the thermodynamic data, typically output from 
            ConcentrationAnalysis and PressureAnalysis class.
            Keys of the dict are set as attributes of ThermoData.
            The items of the dict usually can be:
                partial_pressures : (list)
                    List of partial pressure values.
                variable_defect_specie : (str)
                    Name of variable defect species.
                variable_concentrations : (list)
                    List of concentrations of variable species. 
                defect_concentrations : (list)
                    List of DefectConcentrations objects.
                carrier_concentrations : (list)
                    List of tuples with intrinsic carriers concentrations (holes,electrons).
                conductivities : (list)
                    List of conductivity values (in S/m).
                fermi_levels : (list)
                    list of Fermi level values.
        temperature : (float), optional
            Temperature at which the data is computed. The default is None.
        name : (str), optional
            Name to assign to ThermoData. The default is None.
        """
        self.data = thermodata  
        self.temperature = temperature if temperature else None
        self.name = name if name else None
        for k,v in thermodata.items():
            setattr(self, k, v)
        
    
    def __getitem__(self,key):
        return self.data[key]
    
    def __iter__(self):
        return self.data.__iter__()
    
    def keys(self):
        return self.data.keys()
    
    def values(self):
        return self.data.values()
    
    def items(self):
        return self.data.items()
    
    
    def as_dict(self):
        """
        Returns:
            Json-serializable dict representation of ThermoData
        """
        d = {"@module": self.__class__.__module__,
              "@class": self.__class__.__name__,
              "thermodata":self.data.copy(),    
              "temperature": self.temperature,
              "name": self.name
              }
        if 'defect_concentrations' in self.data.keys():
            d['thermodata']['defect_concentrations'] = [dc.as_dict() for dc in self.defect_concentrations]
        d = jsanitize(d) #convert numpy.float64 to float
        return d        

    def to_json(self,path=None):
        """
        Save ThermoData object as json string or file

        Parameters
        ----------
        path : (str), optional
            Path to the destination file.  If None the name of ThermoData is used as filename.

        Returns
        -------
        d : (str)
            If path is not set a string is returned.
        """
        if not path:
            path = op.join(os.getcwd(),f'thermodata_{self.name}.json')
        d = self.as_dict()
        with open(path,'w') as file:
            json.dump(d,file)
        return
 

    @classmethod
    def from_dict(cls,d):
        """
        Reconstitute a ThermoData object from a dict representation created using
        as_dict().

        Args:
            d (dict): dict representation of ThermoData.

        Returns:
            ThermoData object
        """
        if 'thermodata' in d.keys():
            data = d['thermodata'].copy()
            if 'defect_concentrations' in data.keys():
                data['defect_concentrations'] = [DefectConcentrations.from_dict(dc) for dc in data['defect_concentrations']]
            temperature = d['temperature']
            name = d['name']
        else:
            data = d # recover old 
            temperature = None
            name = None
        return cls(data,temperature,name)
     
         
    @staticmethod
    def from_json(path_or_string):
        """
        Build ThermoData object from json file or string.

        Parameters
        ----------
        path_or_string : (str)
            If an existing path to a file is given the object is constructed reading the json file.
            Otherwise it will be read as a string.

        Returns
        -------
        ThermoData object.

        """
        if op.isfile(path_or_string):
            with open(path_or_string) as file:
                d = json.load(file)
        else:
            d = json.loads(path_or_string)
        return ThermoData.from_dict(d)
            
    
    def get_specific_pressures(self,p_values):
        """
        Get ThermoData object only for specific pressure values. The closest pressure value
        present in list is chosen for each provided value.

        Parameters
        ----------
        p_values : (list)
            List of partial pressure values.

        Returns
        -------
        ThermoData object
        """
        seldata = {}
        for p in p_values:
            pressures = self.partial_pressures
            p_sel = min(pressures, key=lambda x:abs(x-p))
            index = pressures.index(p_sel)
            for k,v in self.data.items():
                for e in v:
                    if v.index(e) == index:
                        if k not in seldata.keys():
                            seldata[k] = []
                        seldata[k].append(e)
        
        temperature = self.temperature if self.temperature else None
        

        name = self.name + 'p_' + '-'.join([str(p) for p in p_values]) if self.name else None
            
        return ThermoData(seldata,temperature=temperature,name=name)
        
    
    def set_data(self,key,value):
        """
        Set data dictionary (chosen over __setitem__ to not override data by accident).
        """
        self.data[key] = value
        setattr(self, key, value)
        return

