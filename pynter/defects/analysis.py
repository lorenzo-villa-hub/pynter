# module for analysing defect calculations


import numpy as np
from scipy.optimize import bisect
from pymatgen.electronic_structure.dos import FermiDos
import matplotlib
import matplotlib.pyplot as plt
from pynter.defects.pmg.pmg_dos import FermiDosCarriersInfo
from pynter.tools.utils import get_object_feature
from pynter.defects.defects import DefectName
from monty.json import MontyDecoder, MSONable
import pandas as pd
import os.path as op
import json


class DefectsAnalysis:
    
    def __init__(self, entries, vbm, band_gap, sort_entries = True, occupation_function='MB'):
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
            occupation_function : (str) 
                Function to compute defect occupation. 
                "FD" for Fermi-Dirac, "MB" for Maxwell-Boltzmann.
        """
        self.entries = self.sort_entries(entries) if sort_entries else entries
        self.vbm = vbm
        self.band_gap = band_gap
        self.occupation_function = occupation_function
        self.groups = self._group_entries()
        self.names = list(self.groups.keys())
        

    def __str__(self):     
        return self.get_dataframe().__str__()
            
    def __repr__(self):
        return self.__str__()
    
    def __iter__(self):
        return self.entries.__iter__()
    
    def _group_entries(self):
        groups = {}
        for e in self.entries:
            if e.name not in groups.keys():
                groups[e.name] = [e]
            else:
                groups[e.name].append(e)
        return groups
            
    
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
            d = json.loads(path_or_string)
        return DefectsAnalysis.from_dict(d)


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
            dtype = sd.defect_type
            dspecie = sd.name.dspecie
            dname = self.select_entries(types=[dtype],defect_specie=dspecie)[0].name
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
        
        charge_transition_levels = {}
        # starting point is first E value
        previous_stable_charges = self.stable_charges(None,fermi_level=energy_range[0],entries=entries)
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
            Dictionary of chemical potentials ({Element: number})   
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
                c = e.defect_concentration(self.vbm, chemical_potentials,temperature,fermi_level,per_unit_volume,self.occupation_function)
                corr = self._get_frozen_correction(e,frozen,dc)
                c = c * corr
                defconc = SingleDefConc(name=e.name,charge=e.charge,conc=c,stable=bool(c<=nsites))
                concentrations.append(defconc)     
            
            else:
                c = e.defect_concentration(self.vbm, chemical_potentials,temperature,fermi_level,per_unit_volume,self.occupation_function)
                defconc = SingleDefConc(name=e.name,charge=e.charge,conc=c,stable=bool(c<=nsites))
                concentrations.append(defconc)
            

        return DefectConcentrations(concentrations)


    def _get_frozen_correction(self,e,frozen,dc):
        corr = 1
        for dn in e.name:
            typ, specie, name = dn.dtype, dn.dspecie, dn.name
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

    
    def filter_entries(self,inplace=False,exclude=False,mode='and',entries=None,
                       types=None,elements=None,names=None,**kwargs):
        """
        Filter entries based on different criteria. Return another DefectsAnalysis object.

        Parameters
        ----------
        inplace : (bool)
            Apply changes to current DefectsAnalysis object.
        exclude : (bool), optional
            Exclude the entries satisfying the criteria instead of selecting them. The default is False.
        mode : (str), optional
            Filtering mode, possibilities are: 'and' and 'or'. The default is 'and'.
        entries : (list)
            List of defect entries.
        types : (list), optional
            Class name of the defect in the entry.
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
        output_entries = self.select_entries(exclude=exclude,mode=mode,entries=entries,types=types,
                                             elements=elements,names=names,**kwargs)
        
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
        

    def get_dataframe(self,entries=None,pretty=False,include_bulk=False,display=[]):
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
    
                
    def plot(self,chemical_potentials,entries=None,xlim=None,ylim=None,title=None,
             fermi_level=None,plotsize=1,fontsize=1.2,show_legend=True,
             format_legend=True,get_subplot=False,subplot_settings=None):
        """
        Produce defect Formation energy vs Fermi energy plot
        Args:
            chemical_potentials:
                Chempots object containing the chemical
                potential of each element.
            entries: 
                List of entries to calculate. If None all entries are considered.
            xlim:
                Tuple (min,max) giving the range of the x (fermi energy) axis.
            ylim:
                Tuple (min,max) giving the range for the formation energy axis.
            fermi_level:
                float to plot Fermi energy position.
            plotsize:
                float or tuple. Can be float or tuple for different x and y sizes
                multiplier to change plotsize.
            fontsize:
                float  multiplier to change fontsize.
            show_legend:
                Bool for showing legend.
            format_legend:
                Bool for getting latex-like legend based on the name of defect entries.
            get_subplot:
                Bool for getting subplot.
            subplot_settings:
                List with integers for subplot setting on matplotlib (plt.subplot(nrows,ncolumns,index)). 
        Returns:
            matplotlib object
        """
        entries = entries if entries else self.entries
        matplotlib.rcParams.update({'font.size': 10*fontsize})        
        formation_energies = self.formation_energies(chemical_potentials,fermi_level=0,entries=entries)
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
            
        for name in formation_energies:
            energy = np.zeros(len(x))
            emin = np.zeros(len(x))
            x_star = []
            y_star = []
            q_previous = None
            # this is faster than using self.stable_charges
            for i in range(0,npoints):
                emin[i] = 1e40
                for d in formation_energies[name]:
                    q = d[0]
                    e0 = d[1]
                    energy[i] = e0 + q*x[i]
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

            if format_legend:
                label_txt = name.symbol
            else:
                label_txt = name            

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
     
        
    def plot_binding_energies(self, names=None, xlim=None, ylim=None, size=1,format_legend=True):
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
                if e.defect_type == 'DefectComplex':
                    if e.name not in names:
                        names.append(e.name)        
        # getting binding energy at different fermi levels for every name in list
        for name in names:
            label = self.select_entries(names=[name])[0].symbol if format_legend else name
            for i in range(0,len(ef)):
                binding_energy[i] = self.binding_energy(name,fermi_level=ef[i])            
            plt.plot(ef,binding_energy, linewidth=2.5*size,label=label)
            
        plt.axvline(x=0.0, linestyle='-', color='k', linewidth=2)  # black lines for gap edges
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
    
    
    def plot_ctl(self, entries=None,ylim=None, size=1, fermi_level=None, format_legend=True,get_integers=True):
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
        plt.figure(figsize=(10*size,10*size))         
        if ylim == None:
            ylim = (-0.5,self.band_gap +0.5)        
        charge_transition_levels = self.charge_transition_levels(entries=entries,get_integers=get_integers)
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
                x_ticks_labels[x_ticks_labels.index(name)] = name.symbol               
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
    
    
    def select_entries(self,exclude=False,mode='and',entries=None,
                       types=None,elements=None,names=None,**kwargs):
        """
        Find entries based on different criteria. Returns a list of DefectEntry objects.

        Parameters
        ----------
        exclude : (bool), optional
            Exclude the entries satisfying the criteria instead of selecting them. The default is False.
        mode : (str), optional
            Filtering mode, possibilities are: 'and' and 'or'. The default is 'and'. 
        entries : (list)
            List of defect entries.
        types : (list), optional
            Class name of the defect in the entry.
        elements : (list), optional
            List of symbols of elements that need to belong to the defect specie.
            If None this criterion is ignored. The default is None.
        names : (list)
            List of entry names.
        **kwargs : (dict)
            Criterion based on attributes or methods of defect entry (e.g. charge=1).

        Returns
        -------
        List of DefectEntry objects.
        """        
        input_entries = entries if entries else self.entries 
        ent = input_entries.copy()
        sel_entries = []
        if types:
            if sel_entries and mode=='and':
                ent = sel_entries.copy()
                sel_entries = []
            for e in ent:
                if e.defect_type in types:
                    sel_entries.append(e)
        
        if elements is not None:
            if sel_entries and mode=='and':
                ent = sel_entries.copy()
                sel_entries = []
            for e in ent:
                select = False                
                for dn in e.name:
                    if dn.dspecie in elements:
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
    
        return output_entries    
    

    def solve_fermi_level(self,chemical_potentials,bulk_dos,temperature=300,
                                    fixed_concentrations=None,external_defects=[],xtol=1e-05):
        """
        Solve charge neutrality and get the value of Fermi level at thermodynamic equilibrium.
        
        Parameters
        ----------
        chemical_potentials : (Dict)
            Dictionary of chemical potentials in the format {Element('el'):chempot}.
        bulk_dos : (CompleteDos object)
            Pymatgen CompleteDos object of the DOS of the bulk system.
        temperature : (float), optional
            Temperature to equilibrate the system to. The default is 300.
        fixed_concentrations: (dict)
            Dictionary with fixed concentrations. Keys are defect entry names in the standard
            format, values are the concentrations (ex {'Vac_A':1e20}) . For more info, 
            read the documentation of the defect_concentrations method.
        external_defects : (list)
            List of SingleDefConc objects containing external defect concentrations 
            (not present in defect entries).
        xtol: Tolerance for bisect (scipy) to solve charge neutrality. The default is 1e-05.

        Returns
        -------
        (float)
            Fermi level dictated by charge neutrality.
        """
        
        fdos = FermiDosCarriersInfo(bulk_dos, bandgap=self.band_gap)
        _,fdos_vbm = fdos.get_cbm_vbm()

        def _get_total_q(ef):
            qd_tot = sum([
                d.charge * d.conc
                for d in self.defect_concentrations(chemical_potentials,temperature,ef,
                                                    fixed_concentrations)])
            for d_ext in external_defects:
                qd_tot += d_ext.charge * d_ext.conc                
            qd_tot += fdos.get_doping(fermi_level=ef + fdos_vbm, temperature=temperature)
            return qd_tot
                       
        return bisect(_get_total_q, -1., self.band_gap + 1.,xtol=xtol)
    
    
    def sort_entries(self,entries=None,features=['name','charge'],inplace=False):
        """
        Sort defect entries with different criteria.

        Parameters
        ----------
        entries : (list), optional
            List of defect entries to sort. If None self.entries is used. The default is None.
        features : (list), optional
            List of strings with attribute/method names. The default is ['name','charge'].
        inplace : (bool), optional
            Reset the self.entries attibute with sorted entries. Only works if entries input
            is not given. The default is False.

        Returns
        -------
        (list)
            List of defect entries.
        """
        inplace = False if entries else inplace
        entries = entries if entries else self.entries
        def criteria(entry):
            return [get_object_feature(entry,feature) for feature in features]           
        sort_function = lambda entry: criteria(entry)
        sorted_entries = sorted(entries, key=sort_function)
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
        

    def _get_total_charge(self,fermi_level,chemical_potentials, bulk_dos, temperature=300, 
                          fixed_concentrations=None,external_defects=[], per_unit_volume=True):
        fdos = FermiDos(bulk_dos, bandgap=self.band_gap)
        _,fdos_vbm = fdos.get_cbm_vbm()    
        qd_tot = sum([
            d.charge * d.conc
            for d in self.defect_concentrations(chemical_potentials,temperature,fermi_level,
                                                fixed_concentrations,per_unit_volume)])
        for d_ext in external_defects:
            qd_tot += d_ext.charge * d_ext.conc
        qd_tot += fdos.get_doping(fermi_level=fermi_level + fdos_vbm, temperature=temperature)
        return qd_tot



class SingleDefConc(MSONable):
    
    def __init__(self,name,charge,conc,stable=True):
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
        stable : (bool)
            If formation energy is positive, i.e. the concentration < Nsites.
        """
        self.name = name
        self.charge = charge
        self.conc = conc
        self.stable = stable
        
    def __repr__(self):
        s = 'charge=%.1f, conc=%.2e, name=%s, stable=%s' %(self.charge,self.conc,self.name,self.stable)
        return s
    
    def __print__(self):
        return self.__repr__()
    
    def __getitem__(self,key):
        return getattr(self,key)
        
        
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
            cmax = SingleDefConc(name='',conc=None,charge=0,stable=True) # dummy object
            for c in concs:
                if cmax.conc is None or c.conc >= cmax.conc:
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
        dc = [SingleDefConc.from_dict(c) for c in d]
        return cls(dc)
 
    @property
    def elemental(self):
        """
        Dictionary with element (or element vacancy) as keys and total element concentration as values.
        """
        d = {}
        for c in self:
            for dn in c.name:
                if dn.dspecie not in d.keys():
                    if dn.dtype == 'Vacancy':
                        ekey = dn.name
                        d[ekey] = self.get_element_total(dn.dspecie,vacancy=True)
                    else:
                        ekey = dn.dspecie
                        d[ekey] = self.get_element_total(dn.dspecie,vacancy=False)
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
                              charges=None,indexes=None,**kwargs):
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
        **kwargs : (dict)
            Criteria for selection. They need to be attributes of SingleDefConc.

        Returns
        -------
        output_concs : (list)
            Filtered DefectConcentrations object.
        """

        output_concs = self.select_concentrations(mode=mode,exclude=exclude,names=names,
                                                  charges=charges,indexes=indexes,**kwargs)
                    
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
            for dn in c.name:
                if element == dn.dspecie:
                    if vacancy and dn.dtype=='Vacancy':
                        eltot += c.conc
                    elif vacancy == False:
                        eltot += c.conc
        return eltot
        

    def select_concentrations(self,concentrations=None,mode='and',exclude=False,names=None,
                    charges=None,indexes=None,**kwargs):
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
        **kwargs : (dict)
            Criteria for selection. They need to be attributes of SingleDefConc.

        Returns
        -------
        output_concs : (list)
            List of SingleDefConc objects.
        """
        sel_concs = []
        input_concs = concentrations if concentrations else self.concentrations.copy()
        concs = input_concs.copy()  
        
        if names is not None:
            if sel_concs and mode=='and':
                concs = sel_concs.copy()
                sel_concs = []
            for c in concs:
                if c.name in names:
                    if c not in sel_concs:
                        sel_concs.append(c)
                        
        if charges is not None:
            if sel_concs and mode=='and':
                concs = sel_concs.copy()
                sel_concs = []
            for c in concs:
                if c.charge in charges:
                    if c not in sel_concs:
                        sel_concs.append(c)
                        
        if indexes is not None:
            if sel_concs and mode=='and':
                concs = sel_concs.copy()
                sel_concs = []
            for c in concs:
                if concs.index(c) in indexes:
                    if c not in sel_concs:
                        sel_concs.append(c)
        
        if kwargs:
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
                    
        
        
        
