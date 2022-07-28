#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:48:55 2021

@author: villa
"""

from pynter.defects.analysis import DefectsAnalysis, DefectConcentrations
from pynter.phase_diagram.experimental import ChempotExperimental
from pynter.tools.utils import save_object_as_json, get_object_from_json
import copy
import os.path as op
import json
from monty.json import MontyDecoder
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
            Dictionary with mobility values for the defect species. 
            Keys must contain "electrons", "holes" and the defect specie name.
        """
        self.mobilities = mobilities
        
        
    def get_conductivity(self,carrier_concentrations,defect_concentrations,temperature=300):
        """
        Calculate conductivity from the concentrations of electrons, holes and defects and their mobilities.
        
        Parameters
        ----------
        carrier_concentrations : (list)
            List of tuples with intrinsic carriers concentrations (holes,electrons).
        defect_concentrations : (list)
            Defect concentrations in the same format as the output of DefectsAnalysis. 
        temperature : float, optional
            Value of temperature. The default is 300.

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
            sigma_ionic += mob[dname] * d['conc'] * abs(d['charge']) * e *1e06 #concentrations need to be in 1/m**3
        sigma = sigma_el + sigma_ionic
        
        return sigma


class PressureAnalysis:
    """
    Class that handles the analysis of the oxygen partial pressure dependency.
    """
    
    def __init__(self,defects_analysis,bulk_dos,frozen_defect_concentrations=None,external_defects=[],xtol=1e-05):
        """
        Parameters
        ----------
        defects_analysis :
            DefectsAnalysis object.
        bulk_dos : 
            Pymatgen Dos object.
        frozen_defect_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations. (ex {'Vac_Na':1e20}) 
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
        xtol: Tolerance for bisect (scipy) to solve charge neutrality. The default is 1e-05.
        """
        self.da = defects_analysis
        self.bulk_dos = bulk_dos
        self.frozen_defect_concentrations = frozen_defect_concentrations if frozen_defect_concentrations else None
        self.external_defects = external_defects if external_defects else []
        self.xtol = xtol
    
    
    def get_concentrations(self,reservoirs,temperature=None,name=None):
        """
        Calculate defect and carrier concentrations at different oxygen partial pressure values

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
        Thermodata : (thermoData)
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
        dos = self.bulk_dos
        frozen_df = self.frozen_defect_concentrations
        ext_df = self.external_defects
        for r,mu in res.items():
            if frozen_df or ext_df:
                mue = self.da.non_equilibrium_fermi_level(frozen_df,mu,dos,ext_df,temperature=T,xtol=self.xtol)
            else:
                mue = self.da.equilibrium_fermi_level(mu,dos,temperature=T,xtol=self.xtol)
            conc = self.da.defect_concentrations(mu,T,mue,frozen_df)
            carriers = self.da.carrier_concentrations(dos,temperature=T,fermi_level=mue)
            defect_concentrations.append(conc)
            carrier_concentrations.append(carriers)
            fermi_levels.append(mue)
        
        thermodata = {}
        thermodata['partial_pressures'] = partial_pressures 
        thermodata['defect_concentrations'] = defect_concentrations 
        thermodata['carrier_concentrations'] = carrier_concentrations
        thermodata['fermi_levels'] = fermi_levels
        
        thermodata = ThermoData(thermodata,temperature=temperature,name=name)
        
        return thermodata
    
    
    def get_conductivities(self,reservoirs,mobilities,temperature=None,name=None):
        """
        Calculate conductivity as a function of oxygen partial pressure.

        Parameters
        ----------
        reservoirs : (dict, Reservoirs or PressureReservoirs)
            Object with partial pressure values as keys and chempots dictionary as values.
        mobilities : (dict)
            Dictionary with mobility values for the defect species. 
            Keys must contain "electrons", "holes" and the defect specie name.
        temperature : (float), optional
            Temperature in Kelvin. If None self.temperature is used. The default is None.
        name : (str), optional
            Name to assign to ThermoData.

        Returns
        -------
        Thermodata : (thermoData)
            ThermoData object that contains the thermodynamic data:
                partial_pressures : (list)
                    List of partial pressure values.
                conductivities : (list)
                    List of conductivity values (in S/m).
        """      
        res = reservoirs
        if temperature:
            T = temperature
        elif hasattr(res,'temperature'):
            T = res.temperature
        else:
            raise ValueError('Temperature needs to be provided or to be present ad attribute in PressureReservoirs object')
        cnd = Conductivity(mobilities)
        partial_pressures = list(res.keys())
        conductivities = []
        dos = self.bulk_dos
        frozen_df = self.frozen_defect_concentrations
        ext_df = self.external_defects
        for r,mu in res.items():
            if frozen_df or ext_df:
                mue = self.da.non_equilibrium_fermi_level(frozen_df,mu,dos,ext_df,temperature=T,xtol=self.xtol)
            else:
                mue = self.da.equilibrium_fermi_level(mu,dos,temperature=T,xtol=self.xtol)
            conc = self.da.defect_concentrations(mu,T,mue,frozen_df)
            carriers = self.da.carrier_concentrations(dos,temperature=T,fermi_level=mue)
            sigma = cnd.get_conductivity(carriers, conc)
            conductivities.append(sigma)
        
        thermodata={}
        thermodata['partial_pressures'] = partial_pressures
        thermodata['conductivities'] = conductivities
        
        thermodata = ThermoData(thermodata,temperature=temperature,name=name)
        
        return thermodata
    
    
    def get_fermi_levels(self,reservoirs,temperature=None,name=None):
        """
        Calculate defect and carrier concentrations at different oxygen partial pressure values

        Parameters
        ----------
        reservoirs : (dict, Reservoirs or PressureReservoirs)
            Object with partial pressure values as keys and chempots dictionary as values.
        temperature : (float), optional
            Temperature in Kelvin. If None self.temperature is used. The default is None.
        name : (str), optional
            Name to assign to ThermoData.
            
        Returns
        -------
        Thermodata : (ThermoData)
            ThermoData object that contains the thermodynamic data:
                partial_pressures : (list)
                    List of partial pressure values.
                fermi_levels : (list)
                    List of Fermi level values
        """
        res = reservoirs
        if temperature:
            T = temperature
        elif hasattr(res,'temperature'):
            T = res.temperature
        else:
            raise ValueError('Temperature needs to be provided or to be present ad attribute in PressureReservoirs object')
        partial_pressures = list(res.keys())
        fermi_levels = []
        dos = self.bulk_dos
        frozen_df = self.frozen_defect_concentrations
        ext_df = self.external_defects
        for r,mu in res.items():
            if frozen_df or ext_df:
                mue = self.da.non_equilibrium_fermi_level(frozen_df,mu,dos,ext_df,temperature=T,xtol=self.xtol)
            else:
                mue = self.da.equilibrium_fermi_level(mu,dos,temperature=T,xtol=self.xtol)
            fermi_levels.append(mue)
        
        thermodata = {}
        thermodata['partial_pressures'] = partial_pressures
        thermodata['fermi_levels'] = fermi_levels
        
        thermodata = ThermoData(thermodata,temperature=temperature,name=name)
        
        return thermodata
    

    def get_quenched_fermi_levels(self,reservoirs,initial_temperature,final_temperature,
                                  quenched_species=None,get_final_concentrations=True,
                                  quench_elements=False,name=None):
        """
        Calculate Fermi level as a function of oxygen partial pressure with quenched defects. 
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
        name : (str), optional
            Name to assign to ThermoData.
        get_final_concentrations : (bool)
            Save also defect and carrier concentrations after quenching. The default is True.
        quench_elements : (bool)
            If True the total concentrations of elements at high temperature go in the charge neutrality at low temperature.
            If False the quenched concentrations are the ones of single defect species (e.g. elements are not allowed
            to equilibrate on different sites). The default is False.

        Returns
        -------
        Thermodata : (thermoData)
            ThermoData object that contains the thermodynamic data:
                partial_pressures : (list)
                    List of partial pressure values.
                fermi_levels : (list)
                    List of Fermi level values
        """
        
        T1 = initial_temperature
        T2 = final_temperature
        res = reservoirs
        if hasattr(res,'temperature'):
            if res.temperature != T1:
                print('Warning: PressureReservoirs temperature is not set to the initial quenching temperature')
        partial_pressures = list(res.keys())
        fermi_levels = []
        if get_final_concentrations:
            final_defect_conc = []
            final_carrier_conc = []
            
        dos = self.bulk_dos
        frozen_df = self.frozen_defect_concentrations
        ext_df = self.external_defects
                
        for r,mu in res.items():
            if frozen_df or ext_df:
                mue = self.da.non_equilibrium_fermi_level(frozen_df,mu,dos,ext_df,temperature=T1,xtol=self.xtol)
            else:
                mue = self.da.equilibrium_fermi_level(mu,dos,temperature=T1,xtol=self.xtol)
            if quench_elements:
                c1 = self.da.defect_concentrations(mu,T1,mue,frozen_df).elemental
            else:
                c1 = self.da.defect_concentrations(mu,T1,mue,frozen_df).total
            if quenched_species is None:
                quenched_concentrations = c1.copy()
            else:
                quenched_concentrations = copy.deepcopy(frozen_df) if frozen_df else {}
                for k in quenched_species:
                    quenched_concentrations[k] = c1[k]
            quenched_mue = self.da.non_equilibrium_fermi_level(quenched_concentrations,mu,dos,ext_df,temperature=T2,xtol=self.xtol)
            fermi_levels.append(quenched_mue)
            
            conc = self.da.defect_concentrations(mu,T2,quenched_mue,quenched_concentrations)
            carriers = self.da.carrier_concentrations(dos,temperature=T2,fermi_level=quenched_mue)
            final_defect_conc.append(conc)
            final_carrier_conc.append(carriers)
                
            
        thermodata = {}
        thermodata['partial_pressures'] = partial_pressures
        thermodata['fermi_levels'] = fermi_levels
        if get_final_concentrations:
            thermodata['defect_concentrations'] = final_defect_conc
            thermodata['carrier_concentrations'] = final_carrier_conc
        
        thermodata = ThermoData(thermodata,temperature=(initial_temperature,final_temperature),name=name)
        
        return thermodata


class ThermoData:
    
    "Class to handle defect thermodynamics data"
    
    def __init__(self,thermodata,temperature=None,name=None):
        """
        Class that handles dict of defect thermodynamics data.
        
        Parameters
        ----------
        thermodata : (dict)
            Dict that contains the thermodynamic data, typically output from PressureAnalysis class.
            Keys of the dict are set as attributes of ThermoData.
            The items of the dict usually can be:
                partial_pressures : (list)
                    List of partial pressure values.
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
        d['thermodata']['defect_concentrations'] = [dc.as_dict() for dc in self.defect_concentrations]
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
            d = json.load(path_or_string)
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
        