#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 16:33:50 2025

@author: lorenzo
"""

import warnings
import json
import os.path as op
from pandas import DataFrame
from monty.json import MSONable, MontyEncoder
from pynter.tools.format import format_composition
import copy

from pymatgen.analysis.phase_diagram import PhaseDiagram

from .core import Chempots
from .phase_diagram import PDHandler


class Reservoirs(MSONable):
    
    def __init__(self,res_dict,phase_diagram=None,mu_refs=None,are_chempots_delta=False):
        """
        Class to handle dictionaries of chemical potentials. Works almost completely like a python dictionary.

        Parameters
        ----------
        res_dict : (dict)
            Dictionary with reservoir names as key and dictionaries of chemical potentials as values.
        phase_diagram : (PhaseDiagram object), optional
            PhaseDiagram object (Pymatgen), useful to convert absolute chempots in referenced chempots. The default is None.
        mu_refs : (dict)
            Dictionary with chemical potentials of reference elements ({Element:chempot}). If the PhaseDiagram is provided
            mu_refs is taken from the mu_refs attribute. 
        are_chempots_delta : (bool), optional
            Set this variable to True if chempots in dictionary are referenced values. The default is False.
        """
        self.res_dict = res_dict
        self.pd = phase_diagram 
        if mu_refs:
            self.mu_refs = mu_refs
        elif self.pd:
            self.mu_refs = PDHandler(self.pd).get_chempots_reference()
        else:
            self.mu_refs = mu_refs
            warnings.warn('Neither PhaseDiagram or reference chempots have been provided, conversions btw ref and abs value will not be possible',UserWarning)
        self._are_chempots_delta = are_chempots_delta


    def __str__(self):
        df = self.get_dataframe()
        return df.__str__()
    
    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self.res_dict)

    def __iter__(self):
        return self.res_dict.keys().__iter__()

    def __getitem__(self,reskey):
        return self.res_dict[reskey]

    def __setitem__(self,reskey,chempots):
        self.res_dict[reskey] = chempots
        return
    
    def __eq__(self, other):
        if isinstance(other, dict):
            return self.res_dict == other
        elif isinstance(other, Reservoirs):
            return self.res_dict == other.res_dict
        else:
            return False

    def keys(self):
        return self.res_dict.keys()

    def values(self):
        return self.res_dict.values()
    
    def items(self):
        return self.res_dict.items()

    def copy(self):
        return Reservoirs(copy.deepcopy(self.res_dict),phase_diagram=self.pd,
                          are_chempots_delta=self.are_chempots_delta,mu_refs=self.mu_refs)
    
    def update(self, other):
        if isinstance(other, dict):
            for key, value in other.items():
                self.res_dict[key] = value
        else:
            for key, value in other:
                self.res_dict[key] = value


    @property
    def are_chempots_delta(self):
        return self._are_chempots_delta

    def as_dict(self):
        """
        Json-serializable dict representation of a Reservoirs object. 
        
        Returns
        -------
        dict
            Json-serializable dict of a Reservoirs object.
        """
        d = {}
        d['@module'] = self.__class__.__module__
        d['@class'] = self.__class__.__name__
        d['res_dict'] = {r:mu.as_dict() for r,mu in self.res_dict.items()} 
        d['phase_diagram'] = self.pd.as_dict() if self.pd else None
        d['mu_refs'] = self.mu_refs.as_dict() if self.mu_refs else None
        d['are_chempots_delta'] = self.are_chempots_delta
        return d


    def to_json(self,path,cls=MontyEncoder):
        """
        Save Reservoirs object as json string or file

        Parameters
        ----------
        path : (str), optional
            Path to the destination file.  If None a string is exported.
        cls : (cls)
            Encoder class for json.dump. The default is MontyEncoder.

        Returns
        -------
        d : (str)
            If path is not set a string is returned.
        """
        d = self.as_dict()
        if path:
            with open(path,'w') as file:
                json.dump(d,file,cls=cls)
            return
        else:
            return d.__str__()  


    @classmethod
    def from_dict(cls,d):
        """
        Constructor of Reservoirs object from dictionary representation.
        
        Parameters
        ----------
        d : dict
        
        Returns
        -------
        Reservoirs object.
        """
        res_dict = {}
        for res,chempots in d['res_dict'].items():
            res_dict[res] = Chempots.from_dict(chempots)
        phase_diagram = PhaseDiagram.from_dict(d['phase_diagram']) if d['phase_diagram'] else None
        mu_refs = Chempots.from_dict(d['mu_refs']) if d['mu_refs'] else None
        are_chempots_delta = d['are_chempots_delta']
            
        return cls(res_dict,phase_diagram,mu_refs,are_chempots_delta)


    @staticmethod
    def from_json(path_or_string):
        """
        Build Reservoirs object from json file or string.

        Parameters
        ----------
        path_or_string : (str)
            If an existing path to a file is given the object is constructed reading the json file.
            Otherwise it will be read as a string.

        Returns
        -------
        Reservoir object.

        """
        if op.isfile(path_or_string):
            with open(path_or_string) as file:
                d = json.load(file)
        else:
            d = json.loads(path_or_string)
        return Reservoirs.from_dict(d)


    def filter_reservoirs(self,inplace=False,elements=None):
        """
        Get new Reservoir object filtering the chempots dictionary.

        Parameters
        ----------
        inplace : (bool)
            Apply changes to current Reservoirs object.
        elements : (list), optional
            List of element symbols. The default is None.

        Returns
        -------
        res : 
            Reservoirs object.
        """
        res = self.copy()
        mu_refs = self.mu_refs.copy()
        filtered_dict = res.res_dict
        if elements:
            d = filtered_dict.copy()
            for r in list(d):
                for el in list(d[r]):
                    if el not in elements:
                        del filtered_dict[r][el]
                        if el in mu_refs.keys():
                            del mu_refs[el]
        if inplace:
            self.res_dict = filtered_dict
            self.mu_refs = mu_refs
            return
        else:
            res.mu_refs = mu_refs
            return res


    def get_absolute_res_dict(self):
        """ 
        Return values of chempots from referenced to absolute 
        """
        if self.are_chempots_delta:
            return {r:mu.get_absolute(self.mu_refs) for r,mu in self.res_dict.items()}
        else:
            raise ValueError('Chemical potential values are already absolute')
                

    def get_referenced_res_dict(self):
        """ 
        Return values of chempots from absolute to referenced 
        """
        if self.are_chempots_delta:
            raise ValueError('Chemical potential values are already with respect to reference')
        else:
            return {r:mu.get_referenced(self.mu_refs) for r,mu in self.res_dict.items()}

    
    def get_dataframe(self,format_symbols=False,format_compositions=False,all_math=False,ndecimals=None):
        """
        Get DataFrame object of the dictionary of reservoirs

        Parameters
        ----------
        format_symbols : (bool), optional
            Format labels of element chempots as \Delta \mu_{\text{"el"}}. The default is False.
        format_compositions : (bool), optional
            Get Latex format of compositions. The default is False.
        all_math : (bool), optional
            Get all characters in composition written in Latex's math format. The default is False.
        ndecimals : (int), optional
            Number of decimals to round the chemical potentials, if None the numbers are not changed.
            The default is None.

        Returns
        -------
        df : 
            DataFrame object.
        """
        res = self._get_res_dict_with_symbols(format_symbols)
        df = DataFrame(res)
        df = df.transpose()
        if format_compositions:
            new_index = []
            for string in df.index:
                new_string = format_composition(string,all_math=all_math)
                new_index.append(new_string)
            df.index = new_index
        if ndecimals:
            df = df.round(decimals=ndecimals)
        return df
        
    
    def get_latex_table(self,ndecimals=1):
        """
        Get string with formatted latex table of chemical potentials
        """
        df = self.get_dataframe(format_symbols=True,ndecimals=ndecimals)
        table = df.to_latex(escape=False)
        return table
    
    
    def get_plot(self,elements,size=1,**kwargs):
        """
        Plot the stability diagram with the reservoir points

        Parameters
        ----------
        elements : (list)
            List of strings with element symbols on the diagram axis.
        size : (float), optional
            Size of the points. The default is 1.
        **kwargs : (dict)
            Kwargs for the add_reservoirs function.

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        from pynter.phase_diagram.plotter import PDPlotterAdder # import here to avoid circular import
        res = self.copy()
        if not res.pd:
            raise ValueError('PhaseDiagram object needs to be stored to plot the stability diagram')
        plt = PDHandler(res.pd).get_stability_diagram(elements)
        if not res.are_chempots_delta:
            res.set_to_referenced()
            
        plt = PDPlotterAdder(res.pd,size).add_reservoirs(res,elements,**kwargs)
        return plt
                
    
    def set_to_absolute(self):
        """
        Set reservoir dictionary to absolute values of chempots.
        """
        new_res_dict = self.get_absolute_res_dict()
        self.res_dict = new_res_dict
        self._are_chempots_delta = False
        return

        
    def set_to_referenced(self):
        """
        Set reservoir dictionary to referenced values of chempots.
        """
        new_res_dict = self.get_referenced_res_dict()
        self.res_dict = new_res_dict
        self._are_chempots_delta = True
        return    
    
                
    def _get_res_dict_with_symbols(self,format_symbols=False):
        """
        format_labels : (bool), optional
            Format labels of element chempots as \Delta \mu_{\text{"el"}}. The default is False.
        """
        new_dict = {}
        for res,chempots in self.res_dict.items():
            new_dict[res] = {}
            chempots = chempots.mu #keep just dict for DataFrame
            for el in chempots:
                if format_symbols:
                    if self.are_chempots_delta:
                        label = '$\Delta \mu_{\text{%s}}$' %el
                    else:
                        label = '$\mu_{\text{%s}}$' %el
                else:
                    label = el
                new_dict[res][label] = chempots[el]
        return new_dict

    
class PressureReservoirs(Reservoirs):
    """
    Subclass of Reservoirs which contains temperature information. Useful for partial pressure analysis. 
    """
    def __init__(self,res_dict,temperature=None,phase_diagram=None,mu_refs=None,are_chempots_delta=False):
        super().__init__(res_dict,phase_diagram,mu_refs,are_chempots_delta)
        self.temperature = temperature
        self.pressures = list(self.res_dict.keys())
        
    
    def __eq__(self, other):
        if isinstance(other, dict):
            return self.res_dict == other
        elif isinstance(other, PressureReservoirs):
            return self.res_dict == other.res_dict
        else:
            return False
        
    
    def as_dict(self):
        """
        Json-serializable dict representation of a PressureReservoirs object. The Pymatgen element 
        is expressed with the symbol of the element.

        Returns
        -------
        dict
            Json-serializable dict of a PressureReservoirs object.
        """
        d = {}
        d['@module'] = self.__class__.__module__
        d['@class'] = self.__class__.__name__
        d['res_dict'] = {r:mu.as_dict() for r,mu in self.res_dict.items()}
        d['temperature'] = self.temperature
        d['phase_diagram'] = self.pd.as_dict() if self.pd else None
        d['mu_refs'] = self.mu_refs.as_dict() if self.mu_refs else None
        d['are_chempots_delta'] = self.are_chempots_delta
        return d


    @classmethod
    def from_dict(cls,d):
        """
        Constructor of PressureReservoirs object from dictionary representation.
        
        Parameters
        ----------
        d : dict
        
        Returns
        -------
        PressureReservoirs object.
        """
        res_dict = {float(r):Chempots.from_dict(mu) for r,mu in d['res_dict'].items()}
        temperature = d['temperature'] if 'temperature' in d.keys() else None
        phase_diagram = PhaseDiagram.from_dict(d['phase_diagram']) if d['phase_diagram'] else None
        mu_refs = Chempots.from_dict(d['mu_refs'])  if d['mu_refs'] else None
        are_chempots_delta = d['are_chempots_delta']
            
        return cls(res_dict,temperature,phase_diagram,mu_refs,are_chempots_delta)


    @staticmethod
    def from_json(path_or_string):
        """
        Build PressureReservoirs object from json file or string.

        Parameters
        ----------
        path_or_string : (str)
            If an existing path to a file is given the object is constructed reading the json file.
            Otherwise it will be read as a string.

        Returns
        -------
        PressureReservoir object.

        """
        if op.isfile(path_or_string):
            with open(path_or_string) as file:
                d = json.load(file)
        else:
            d = json.load(path_or_string)
        return PressureReservoirs.from_dict(d)   
    
    
    