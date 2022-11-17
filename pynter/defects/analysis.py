# module for analysing defect calculations


import numpy as np
from scipy.optimize import bisect
from pymatgen.core.units import kb
from pymatgen.core.structure import Structure, PeriodicSite, Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.dos import FermiDos
import matplotlib
import matplotlib.pyplot as plt
from pynter.defects.pmg_dos import FermiDosCarriersInfo
from pynter.defects.utils import get_delta_atoms
from pynter.defects.entries import SingleDefectEntry, DefectComplexEntry, get_formatted_legend
from pynter.tools.utils import get_object_feature, get_object_from_json
from monty.json import MontyDecoder, MSONable
import pandas as pd
import os.path as op
import json


class DefectsAnalysis:
    """ 
    Class to compute defect properties starting from single calculations of defects
    Args:
        entries (list): A list of SingleDefectEntry or DefectComplexEntry objects.
        vbm (float): Valence Band energy to use for all defect entries.
            NOTE if using band shifting-type correction then this VBM
            should still be that of the GGA calculation
            (the bandedgeshifting_correction accounts for shift's
            contribution to formation energy).
        band_gap (float): Band gap to use for all defect entries.
            NOTE if using band shifting-type correction then this gap
            should still be that of the Hybrid calculation you are shifting to.       
    """    
    def __init__(self, entries, vbm, band_gap, sort_entries = True):
        self.entries = sorted(entries, key=lambda x: x.name) if sort_entries else entries
        self.vbm = vbm
        self.band_gap = band_gap


    def __str__(self):     
        return self.get_dataframe().__str__()
            
    def __repr__(self):
        return self.__str__()
    
    
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
        Reconstitute a DefectsAnalysis object from a dict representation created using
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
    def from_json(path_or_string):
        """
        Build Dataset object from json file or string.

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
            d = json.load(path_or_string)
        return DefectsAnalysis.from_dict(d)
    
    
    @property
    def names(self):
        """
        Returns a list with all the different names of defect entries
        """
        
        names = []
        for d in self.entries:
            if d.name not in names:
                names.append(d.name)
        
        return names
    
    
    def filter_entries(self,exclude=False,mode='and',entries=None,entry_class=None,elements=None,names=None,**kwargs):
        """
        Filter entries based on different criteria. Return another DefectsAnalysis object.

        Parameters
        ----------
        entries : (list)
            List of defect entries.
        entry_class : (str), optional
            Class name of the entry. "SingleDefectEntry" or "DefectComplexEntry".
        elements : (list), optional
            List of symbols of elements that need to belong to the defect specie.
            If None this criterion is ignored. The default is None.
        names : (list)
            List of entry names.
        **kwargs : (dict)
            Criterion based on attributes or methods of defect entry (e.g. charge=1).

        Returns
        -------
        DefectsAnalysis object.

        """
        input_entries = entries if entries else self.entries 
        ent = input_entries.copy()
        sel_entries = []
        if entry_class:
            if sel_entries and mode=='and':
                ent = sel_entries.copy()
                sel_entries = []
            for e in ent:
                if e.classname == entry_class:
                    sel_entries.append(e)
        
        if elements is not None:
            if sel_entries and mode=='and':
                ent = sel_entries.copy()
                sel_entries = []
            for e in ent:
                select = False                
                for d in e.defect_species:
                    if d['specie'] in elements:
                        select = True
                if select:
                    sel_entries.append(e)
        
        if names:
            if sel_entries and mode=='and':
                ent = sel_entries.copy()
                sel_entries = []
            for e in ent:
                if e.name in names:
                    sel_entries.append(e)
        
        for feature in kwargs:
            if sel_entries and mode=='and':
                ent = sel_entries.copy()
                sel_entries = []
            for e in ent:
                attr = get_object_feature(e,feature)
                if attr == kwargs[feature]:
                    sel_entries.append(e)        

        output_entries = []
        for e in input_entries:
            if exclude:
                if e not in sel_entries:
                    output_entries.append(e)
            else:
                if e in sel_entries:
                    output_entries.append(e) 
                    
        return DefectsAnalysis(output_entries, self.vbm, self.band_gap)
    
    
    def find_entries_by_type_and_specie(self,dtype,dspecie):
        """
        Find entries by selecting defect type and defect specie. Finds only SingleDefectEntry objects.

        Parameters
        ----------
        dtype : (str)
            Class name of defect type ('Vacancy','Interstitial','Substitution').
        dspecie : (str)
            Symbol of the element of defect specie.

        Returns
        -------
        sel_entries : (list)
            List of SingleDefectEntry objects.
        """
        sel_entries = []
        for e in self.entries:
            if e.classname == 'SingleDefectEntry':
                if e.defect.__class__.__name__ == dtype and e.defect.site.specie.symbol == dspecie:
                    sel_entries.append(e)
        return sel_entries
    
    def find_entries_by_name(self,name):        
        sel_entries = []
        for e in self.entries:
            if e.name == name:
                sel_entries.append(e)
        return sel_entries
               
    def find_entry_by_name_and_charge(self,name,charge):
        for e in self.find_entries_by_name(name):
            if e.charge == charge:
                return e


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
        entry = self.find_entry_by_name_and_charge(name, charge)
        for sd in entry.defect_list:
            dtype = sd.__class__.__name__
            dspecie = sd.site.specie.symbol
            dname = self.find_entries_by_type_and_specie(dtype,dspecie)[0].name
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
        
        fdos = FermiDosCarriersInfo(bulk_dos, bandgap=self.band_gap)
        _,fdos_vbm = fdos.get_cbm_vbm()    
        
        h, n = fdos.get_doping(fermi_level=fermi_level + fdos_vbm, temperature=temperature,carriers_values=True)
        
        return abs(h) , abs(n)

           
    def charge_transition_level(self,name,q1,q2):
        """
        Args:
            chemical_potentials: 
                a dictionnary of {Element:value} giving the chemical
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
    
    
    def charge_transition_levels(self, energy_range=None):
        """
        Computes charge transition levels for all defect entries
        Args:
            energy_range (Tuple): Energy range to evaluate the charge transition levels in
                                  default to (-0.5, Eg + 0.5)                  
        Returns:
            Dictionary with defect name and list of tuples for charge transition levels:
                {name:[(q1,q2,ctl),(q2,q3,ctl),...}
        """                
        
        # set charge_transition_levels dict
        charge_transition_levels = {}
        for name in self.names:
            charge_transition_levels[name] = []
        if energy_range == None:
            energy_range = (-0.5,self.band_gap +0.5)
        
        # creating energy array
        npoints = 3000
        step = abs(energy_range[1]-energy_range[0])/npoints
        e = np.arange(energy_range[0],energy_range[1],step)
        
        # starting point is first E value
        previous_stable_charges = self.stable_charges(None,fermi_level=energy_range[0])
        for i in range(0,len(e)):
            stable_charges = self.stable_charges(None, fermi_level = e[i])
            for name in stable_charges:
                new_charge = stable_charges[name][0]
                previous_charge = previous_stable_charges[name][0]
                if new_charge != previous_charge:
                    previous_stable_charges = stable_charges
                    charge_transition_levels[name].append((previous_charge,new_charge,e[i]))
                            
        return charge_transition_levels
    

    def _get_frozen_correction(self,e,frozen,dc):
        corr = 1
        for ds in e.defect_species:
            typ, specie, name = ds['type'], ds['specie'], ds['name']
            if typ == 'Vacancy':
                k = f'Vac_{specie}'
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
    

    def defect_concentrations(self, chemical_potentials, temperature=300, fermi_level=0.,
                              frozen_defect_concentrations=None,per_unit_volume=True):
        """
        Give list of all concentrations at specified efermi.
        If frozen_defect_concentration is provided the concentration of defect entries are 
        corrected according to the fixed provided values. reference to paper once (if?) is
        published.
            
        Labels are ignored when accounting for the defect species equilibrium.
            
        Parameters
        ----------
        chemical_potentials: (dict)
            Dictionary of chemical potentials ({Element: number})   
        temperature: (float) 
            Temperature to produce concentrations at.
        fermi_level: (float) 
            Fermi level relative to valence band maximum. The default is 0. 
        frozen_defect_concentrations: (dict)
            Dictionary with fixed concentrations. Keys can be simple element strings
            (or vacancies of elements) in the format 'Vac_{el}') if that element needs to be 
            fixed across all defect species, alternatively defect entry names can be used as well 
            to target specific defect entries. The values are the concentrations.
        per_unit_volume: (bool)
            Get concentrations in cm^-3. If False they are per unit cell. The default is True.
        
        Returns:
        --------
            list of dictionaries of defect concentrations
        """
        concentrations = []
        if frozen_defect_concentrations:
            dc = self.defect_concentrations(chemical_potentials,temperature,fermi_level,
                                            frozen_defect_concentrations=None,per_unit_volume=per_unit_volume)
            frozen = frozen_defect_concentrations 

        for e in self.entries:
            nsites = e.multiplicity * 1e24 / e.bulk_structure.volume if per_unit_volume else e.multiplicity
            # frozen defects approach
            if frozen_defect_concentrations:
                c = e.defect_concentration(self.vbm, chemical_potentials,temperature,fermi_level,per_unit_volume)
                corr = self._get_frozen_correction(e,frozen,dc)
                c = c * corr
                defconc = SingleDefConc(name=e.name,charge=e.charge,conc=c,defect_species=e.defect_species,stable=bool(c<=nsites))
                concentrations.append(defconc)     
            
            else:
                c = e.defect_concentration(self.vbm, chemical_potentials,temperature,fermi_level,per_unit_volume)
                defconc = SingleDefConc(name=e.name,charge=e.charge,conc=c,defect_species=e.defect_species,stable=bool(c<=nsites))
                concentrations.append(defconc)
            

        return DefectConcentrations(concentrations)
            
 
    def equilibrium_fermi_level(self, chemical_potentials, bulk_dos, temperature=300, xtol=1e-05):
        """
        Solve for the Fermi energy self-consistently as a function of T
        Observations are Defect concentrations, electron and hole conc
        Args:
            temperature: Temperature to equilibrate fermi energies for
            chemical_potentials: dict of chemical potentials to use for calculation fermi level
            bulk_dos: bulk system dos (pymatgen Dos object)
            xtol: Tolerance for bisect (scipy) to solve charge neutrality. The default is 1e-05.
        Returns:
            Fermi energy dictated by charge neutrality
        """

        fdos = FermiDos(bulk_dos, bandgap=self.band_gap)
        _,fdos_vbm = fdos.get_cbm_vbm()

        def _get_total_q(ef):

            qd_tot = sum([
                d.charge * d.conc
                for d in self.defect_concentrations(
                    chemical_potentials=chemical_potentials, temperature=temperature, fermi_level=ef)
            ])
            qd_tot += fdos.get_doping(fermi_level=ef + fdos_vbm, temperature=temperature)
            return qd_tot

        return bisect(_get_total_q, -1., self.band_gap + 1.,xtol=xtol)
            
        
    def non_equilibrium_fermi_level(self, frozen_defect_concentrations, chemical_potentials, bulk_dos, 
                                        external_defects=[], temperature=300, xtol=1e-05):
        
        # Since last defect_concentrations update could be integrated easily with equilibrium_fermi_level.
        # I kept 2 different functions to limit confusions and facilitate integrations with old notebooks
        """
        Solve charge neutrality in non-equilibrium conditions (when some concentrations are fixed).
        If frozen_defect_concentration is not None the concentration for every SingleDefectEntry 
        is normalized starting from the input fixed concentration as:
            C = C_eq * (Ctot_fix/Ctot_eq)
        while for DefectComplexEntry this is applied for every single defect specie which
        is composing the complex (needs to be present in computed entries):
            C = C_eq * prod(Ctot_fix(i)/Ctot_eq(i))
            
        If external defects are present their charge concentrations are treated as d['charge']*d['conc'].
        
        Parameters
        ----------
        frozen_defect_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations. (ex {'Vac_Na':1e20}) 
        chemical_potentials : (Dict)
            Dictionary of chemical potentials in the format {Element('el'):chempot}.
        bulk_dos : (CompleteDos object)
            Pymatgen CompleteDos object of the DOS of the bulk system.
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
        temperature : (float), optional
            Temperature to equilibrate the system to. The default is 300.
        xtol: Tolerance for bisect (scipy) to solve charge neutrality. The default is 1e-05.

        Returns
        -------
        (float)
            Fermi level dictated by charge neutrality .
        """
        
        fdos = FermiDosCarriersInfo(bulk_dos, bandgap=self.band_gap)
        _,fdos_vbm = fdos.get_cbm_vbm()

        def _get_total_q(ef):
            
            # get groups D1 and D2
            qd_tot = sum([
                d.charge * d.conc
                for d in self.defect_concentrations(chemical_potentials,temperature,ef,
                                                    frozen_defect_concentrations)])

            #external fixed defects - D3
            for d_ext in external_defects:
                qd_tot += d_ext.charge * d_ext.conc
                
            qd_tot += fdos.get_doping(fermi_level=ef + fdos_vbm, temperature=temperature)
            return qd_tot
                       
        return bisect(_get_total_q, -1., self.band_gap + 1.,xtol=xtol)
            
     
    def get_dataframe(self,filter_names=None,pretty=False,include_bulk=False,display=[]):
        """
        Get DataFrame to display entries. 

        Parameters
        ----------
        filter_names : (list), optional
            Only entries whose name is on the list will be displayed. The default is None.
        include_bulk: (bool), optional
            Include bulk composition and space group for each entry in DataFrame.
        display: (list)
            List of strings with defect entry attributes or method results to display.

        Returns
        -------
        df : 
            DataFrame object.
        """
        if filter_names:
            entries = []
            for name in filter_names:
                en = self.find_entries_by_name(name)
                for e in en:
                    entries.append(e)
        else:
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
    
    
    def formation_energies(self,chemical_potentials,fermi_level=0):
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
        Returns:
            {name: [(charge,formation energy)] }
        """
        computed_charges = {}
        for d in self.entries:
            name = d.name
            if name in computed_charges:
                computed_charges[name].append((d.charge,d.formation_energy(self.vbm,chemical_potentials,fermi_level=fermi_level)))
            else:
                computed_charges[name] = []
                computed_charges[name].append((d.charge,d.formation_energy(self.vbm,chemical_potentials,fermi_level=fermi_level)))
        
        return computed_charges
            
    
    def plot(self,mu_elts=None,filter_names=None,xlim=None, ylim=None, title=None, fermi_level=None, 
             plotsize=1, fontsize=1.2, show_legend=True, format_legend=True, 
             order_legend=False, get_subplot=False, subplot_settings=None):
        """
        Produce defect Formation energy vs Fermi energy plot
        Args:
            mu_elts:
                a dictionnary of {Element:value} giving the chemical
                potential of each element
            filter_names: (list)
                Display plot only of defect specie whose names are in the list.
                If None all defect species are displayed. The default is None.
            xlim:
                Tuple (min,max) giving the range of the x (fermi energy) axis
            ylim:
                Tuple (min,max) giving the range for the formation energy axis
            fermi_level:
                float to plot Fermi energy position
            plotsize:
                float or tuple
                can be float or tuple for different x and y sizes
                multiplier to change plotsize            
            fontsize:
                float  multiplier to change fontsize
            show_legend:
                Bool for showing legend
            format_legend:
                Bool for getting latex-like legend based on the name of defect entries
            order_legend:
                Bool for ordering legend by alphabetical order of names
            get_subplot:
                Bool for getting subplot
            subplot_settings:
                List with integers for subplot setting on matplotlib (plt.subplot(nrows,ncolumns,index)) 
        Returns:
            matplotlib object
        """
        
        matplotlib.rcParams.update({'font.size': 10*fontsize})        
        computed_charges = self.formation_energies(mu_elts,fermi_level=0)
        # sort by alphabetical order of name for the plot
        if order_legend:
            computed_charges = {k: v for k, v in sorted(computed_charges.items(), key=lambda item: item[0])}                
        if xlim == None:
            xlim = (-0.5,self.band_gap+0.5)        
        npoints = 200
        step = abs(xlim[1]+0.1-xlim[0])/npoints
        x = np.arange(xlim[0],xlim[1]+0.1,step)
        try:
            if len(plotsize)==2:
                pass
        except:
            plotsize = (plotsize,plotsize)
        
        if get_subplot:
            if subplot_settings[2] == 1:
                plt.figure(figsize=(8*plotsize[0],6*plotsize[1]))               
            plt.subplot(subplot_settings[0],subplot_settings[1],subplot_settings[2])
            plt.grid()
        else:
            plt.figure(figsize=(8*plotsize[0],6*plotsize[1]))
            plt.grid()
            
        for name in computed_charges:
            energy = np.zeros(len(x))
            emin = np.zeros(len(x))
            x_star = []
            y_star = []
            q_previous = None
            for i in range(0,npoints):
                emin[i] = 1e40
                for d in computed_charges[name]:
                    q = d[0]
                    e0 = d[1]
                    energy[i] = e0 + q*x[i] # Eform(x) = Eform(0) + q*x
                    # finding most stable charge state
                    if energy[i] < emin[i]:
                        emin[i] = energy[i]
                        q_stable = q
                # getting data to plot transition levels        
                if q_stable != q_previous:
                    if q_previous != None:
                        x_star.append(x[i])
                        y_star.append(emin[i])
                    q_previous = q_stable       
            # if format_legend is True get latex-like legend
            if format_legend:
                label_txt = self._get_formatted_legend(name)
            else:
                label_txt = name            
            if filter_names:
                if name in filter_names:
                    plt.plot(x,emin,label=label_txt,linewidth=3)
                    plt.scatter(x_star,y_star,s=120,marker='*')
            else:
                plt.plot(x,emin,label=label_txt,linewidth=3)
                plt.scatter(x_star,y_star,s=120,marker='*')
                        
        plt.axvline(x=0.0, linestyle='-', color='k', linewidth=2)  # black dashed lines for gap edges
        plt.axvline(x=self.band_gap, linestyle='-', color='k',
                    linewidth=2)        
        if fermi_level:
            plt.axvline(x=fermi_level, linestyle='dashed', color='k', linewidth=1.5, label='$\mu _{e}$')                
        # shaded areas
        plt.axvspan(xlim[0], 0, facecolor='k', alpha=0.2)
        plt.axvspan(self.band_gap, xlim[1]+0.1, facecolor='k', alpha=0.2)
        plt.hlines(0,xlim[0],xlim[1]+0.1,colors='k',linestyles='dashed',alpha=0.5)
        plt.xlim(xlim)
        if ylim: 
            plt.ylim(ylim) 
        plt.xlabel('Fermi level (eV)')
        plt.ylabel('Formation energy (eV)')
        if title:
            plt.title(title)
        if show_legend:    
            plt.legend()
     #   plt.grid()

        return plt
     

    def _get_formatted_legend(self,fullname):
            return get_formatted_legend(fullname)
    
    
    def plot_binding_energies(self, names=None, xlim=None, ylim=None, size=1,format_legend=True):
        """
        Plot binding energies for complex of defects as a function of the fermi level
        Args:
            names: 
                List of strings with names of DefectComplexEntry. If None all DefectComplexEntry
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
        plt.figure(figsize=(8*size,6*size))
        matplotlib.rcParams.update({'font.size': 10*1.8*size}) 
        if xlim==None:
            xlim = (-0.5,self.band_gap+0.5)
        # building array for x values (fermi level)    
        ef = np.arange(xlim[0],xlim[1]+0.1,(xlim[1]-xlim[0])/200)        
        binding_energy = np.zeros(len(ef))        
        if not names:
            names = []
            for e in self.entries:
                if e.classname == 'DefectComplexEntry':
                    if e.name not in names:
                        names.append(e.name)        
        # getting binding energy at different fermi levels for every name in list
        for name in names:
            label = self._get_formatted_legend(name) if format_legend else name
            for i in range(0,len(ef)):
                binding_energy[i] = self.binding_energy(name,fermi_level=ef[i])            
            plt.plot(ef,binding_energy, linewidth=2.5*size,label=label)
            
        plt.axvline(x=0.0, linestyle='-', color='k', linewidth=2)  # black dashed lines for gap edges
        plt.axvline(x=self.band_gap, linestyle='-', color='k',
                    linewidth=2)        
        # shaded areas
        plt.axvspan(xlim[0], 0, facecolor='k', alpha=0.2)
        plt.axvspan(self.band_gap, xlim[1], facecolor='k', alpha=0.2)
        plt.hlines(0,xlim[0],xlim[1],colors='k',linestyles='dashed',alpha=0.5)
        plt.legend()
        plt.xlim(xlim)
        if ylim: 
            plt.ylim(ylim) 
        plt.xlabel('Fermi level (eV)')
        plt.ylabel('Binding energy (eV)')
        
        return plt
    
    
    def plot_ctl(self, ylim=None, size=1, fermi_level=None, format_legend=True):
        """
        Plotter for the charge transition levels
        Args:
            ylim (Tuple): y-axis limits
            fermi_level (float) : float to plot Fermi energy position
            size (float) : Float multiplier for plot size
            format_legend (bool): Bool for getting latex-like legend based on the name of defect entries 
        Returns:
            matplotlib object                
        """        
        plt.figure(figsize=(10*size,10*size))         
        if ylim == None:
            ylim = (-0.5,self.band_gap +0.5)        
        charge_transition_levels = self.charge_transition_levels()
        number_defects = len(charge_transition_levels)   
        x_max = 10
        interval = x_max/(number_defects + 1)
        x = np.arange(0,x_max,interval)        
        # position of x labels
        x_ticks_positions = []
        for i in range(0,len(x)-1):
            x_ticks_positions.append((x[i+1]-x[i])/2 + x[i])            
        x_ticks_labels = []
        for name in charge_transition_levels:
            x_ticks_labels.append(name)        
        # draw vertical lines to separte defect types
        for i in x:
            plt.axvline(x=i, linestyle='-', color='k', linewidth=1.2, alpha=1, zorder=1)
        xlim = (x[0],x[-1])        
        #VBM and CBM shaded
        plt.axhspan(ylim[0], 0, facecolor='grey', alpha=0.9, zorder=2)
        plt.axhspan(self.band_gap,ylim[1], facecolor = 'grey', alpha=0.9, zorder=2)                
        # plot CTL
        for i in range(0,len(x_ticks_labels)):
            name = x_ticks_labels[i]
            for ctl in charge_transition_levels[name]:
                energy = ctl[2]
                plt.hlines(energy,x[i],x[i+1],colors='k',linewidth=2.25, zorder=3)
                charge1 = '+' + str(ctl[1]) if ctl[1] > 0 else str(ctl[1])
                charge2 = '+' + str(ctl[0]) if ctl[0] > 0 else str(ctl[0])
                label_charge = '(' + charge2 + '/' + charge1 + ')'
                font_space = abs(ylim[1]-ylim[0]) / 100
                if energy < ylim[1] and energy > ylim[0]:
                    plt.text(x[i]+(interval/2)*2/number_defects ,energy+font_space,label_charge,fontsize=16*size)        
        # format latex-like legend
        if format_legend:    
             for name in x_ticks_labels:            
                x_ticks_labels[x_ticks_labels.index(name)] = self._get_formatted_legend(name)                
        if fermi_level:
            plt.axhline(y=fermi_level, linestyle='dashed', color='k', linewidth=1.5, label='$\mu _{e}$')   
        
        plt.text(x[-1]+interval/8,-0.3,'VB',fontsize=25*size)
        plt.text(x[-1]+interval/8,self.band_gap+0.2,'CB',fontsize=25*size)
        plt.xticks(ticks=x_ticks_positions,labels=x_ticks_labels,fontsize = (25-number_defects)*size)
        plt.tick_params(axis='x',length=0,width=0)
        plt.yticks(fontsize=16*size)
        plt.xlim(xlim)  
        plt.ylim(ylim)
        plt.ylabel('Energy(eV)',fontsize=20*size)  
        
        return plt  
            
    
    def stable_charges(self,chemical_potentials,fermi_level=0):
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
        Returns:
            {name:(stable charge, formation energy)}
       """
        computed_charges = self.formation_energies(chemical_potentials,fermi_level=fermi_level)        
        stable_charges = {}
        for name in computed_charges:
            emin = 1e40
            for d in computed_charges[name]:
                q = d[0]
                energy = d[1]
                # finding most stable charge state
                if energy < emin:
                    emin = energy
                    q_stable = q
            stable_charges[name] = (q_stable,emin)
            
        return stable_charges
        
            

class SingleDefConc:
    
    def __init__(self,**kwargs):
        """
        Object to store defect concentrations data. Is also subscribtable like a dictionary.

        Parameters
        ----------
        name : (str)
            Name of defect specie.
        charge : (int)
            Charge of defect specie.
        conc : (float)
            Concentration value.
        defects_species : (list)
            Dictionary with defect type and defect specie (attribute of defect entry).
        stable : (bool)
            If formation energy is positive, i.e. the concentration < Nsites.
        """
        for k,v in kwargs.items():
            setattr(self, k, v)
        
    def __repr__(self):
        s = ''
        for k,v in sorted(vars(self).items(), key=lambda x: x[0]):
            s += f"{k}={v},  "            
        return s
    
    def __print__(self):
        return self.__repr__()
    
    def __getitem__(self,key):
        return getattr(self,key)
        
    def as_dict(self):
        d = {}
        for k,v in vars(self).items():
            d[k] = v
        return d
    
    @classmethod
    def from_dict(cls,d):
        return cls(**d)

        
class DefectConcentrations:
    
    def __init__(self,concentrations,get_old_format=False):
        """
        Class to store sets of defect concentrations (output of concentration calculations with DefectsAnalysis).
        List of SingleDefConc objects. Subscriptable like a list.

        Parameters
        ----------
        concentrations : (list)
            List of SingleDefConc objects.
        """
        self.concentrations = concentrations
        d = {}
        # calculate total concentrations
        for c in self.concentrations:
            gname = c.name.split('(')[0] #name without label
            if gname not in d.keys():
                d[gname] = 0
            d[gname] += c.conc    
        self._total = d
        # store stable concentrations
        conc_stable = []
        for n in self.names:
            concs = self.filter_conc(name=n)
            cmax = SingleDefConc(name='',conc=0,charge=0) # dummy object
            for c in concs:
                if c.conc > cmax.conc:
                    cmax = c
            conc_stable.append(concs[concs.index(cmax)])
        self._stable = conc_stable
        
        
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

   
    def as_dict(self):
        d = [c.as_dict() for c in self.concentrations]
        return d

    @classmethod 
    def from_dict(cls,d):
        try:
            dc = [SingleDefConc.from_dict(c) for c in d]
        except: # recover old thermodata with just total concentrations
            dc = []
            for k,v in d.items():
                c = SingleDefConc(name=k,charge=100,conc=v) # charge if 100 to make it clear it's a dummy value that can't be used
                dc.append(c)
        return cls(dc)
 
    @property
    def elemental(self):
        """
        Dictionary with element (or element vacancy) as keys and total element concentration as values.
        """
        d = {}
        for c in self:
            for ds in c.defect_species:
                if ds['specie'] not in d.keys():
                    if ds['type'] == 'Vacancy':
                        ekey = 'Vac_' + ds['specie']
                        d[ekey] = self.get_element_total(ds['specie'],vacancy=True)
                    else:
                        ekey = ds['specie']
                        d[ekey] = self.get_element_total(ds['specie'],vacancy=False)
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

    
    def filter_conc(self,mode='and',exclude=False,**kwargs):
        """
        Filter concentrations based on different criteria.

        Parameters
        ----------
        mode : (str), optional
            Selection mode. available are 'and','or'. The default is 'and'.
        exclude : (bool), optional
            Exclude the entries matching the criteria. The default is False.
        **kwargs : (dict)
            Criteria for selection. They need to be attributes of SingleDefConc.

        Returns
        -------
        output_concs : (list)
            Filtered DefectConcentrations object.
        """
        sel_concs = []
        input_concs = self.concentrations.copy()
        concs = input_concs.copy()  
        for k in kwargs:
            if sel_concs and mode=='and':
                concs = sel_concs.copy()
                sel_concs = []
            for c in concs:
                feature = get_object_feature(c,k)
                if feature == kwargs[k]:
                    if c not in sel_concs:
                        sel_concs.append(c)
                        
        output_concs = []
        for c in input_concs:
            if exclude:
                if c not in sel_concs:
                    output_concs.append(c)
            else:
                if c in sel_concs:
                    output_concs.append(c)
                    
        return output_concs
            
        
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
            for ds in c.defect_species:
                if element == ds['specie']:
                    if vacancy and ds['type']=='Vacancy':
                        eltot += c.conc
                    elif vacancy == False:
                        eltot += c.conc
        return eltot
        
        
        
        
        
