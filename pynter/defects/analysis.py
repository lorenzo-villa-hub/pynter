# module for analysing defect calculations


import numpy as np
from scipy.optimize import bisect
from pymatgen.analysis.defects.utils import kb
from pymatgen.core.structure import Structure, PeriodicSite, Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.dos import FermiDos
import matplotlib
import matplotlib.pyplot as plt
from pynter.defects.pmg_dos import FermiDosCarriersInfo
from pynter.defects.utils import get_delta_atoms
from pynter.defects.entries import SingleDefectEntry

    
    
class DefectsAnalysis:
    """ 
    Class to compute defect properties starting from single calculations of defects
    Args:
        entries (list): A list of SingleDefectEntry objects
        vbm (float): Valence Band energy to use for all defect entries.
            NOTE if using band shifting-type correction then this VBM
            should still be that of the GGA calculation
            (the bandedgeshifting_correction accounts for shift's
            contribution to formation energy).
        band_gap (float): Band gap to use for all defect entries.
            NOTE if using band shifting-type correction then this gap
            should still be that of the Hybrid calculation you are shifting to.       
    """    
    def __init__(self, entries, vbm, band_gap):
        self.entries = entries
        self.vbm = vbm
        self.band_gap = band_gap


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
        entries = [SingleDefectEntry.from_dict(entry_dict) for entry_dict in d['entries']]        
        vbm = d['vbm']
        band_gap = d['band_gap']
        return cls(entries,vbm,band_gap)
    
    
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
            name (string): name of defect complex as assigned in SingleDefectData object
            fermi_level (float): Fermi level at which binding energy needs to be computed
        Returns:
            binding_energy (float)
        """
        # finding entry associated to 'name'        
        for entry in self.entries:
            if entry.name == name:
                defect_complex = entry
                print(defect_complex)
                break
        
        # finding names and elements of single defects of which the complex is made
        single_defects = {}
        for el in entry.delta_atoms:
            # sign of delta atoms - if > 0 interstitial, if < 0 vacancy
            sign_el = 1 if entry.delta_atoms[el]>0 else -1 if entry.delta_atoms[el]<0 else 0
            for d in self.entries:
                # has to be a single defect - delta_atoms needs to have only 1 key
                if el in d.delta_atoms and len(d.delta_atoms.keys()) == 1:
                    # sign in delta atoms of single defect has to be the same (Vacancy or interstitial)
                    sign_el_single = 1 if d.delta_atoms[el]>0 else -1 if d.delta_atoms[el]<0 else 0 
                    if sign_el == sign_el_single:
                        # add name to dict of single defects to consider for binding energy
                        single_defects[el] = d.name
        
        # getting stable charge states at desired fermi level                    
        stable_charges = self.stable_charges(None,fermi_level=fermi_level)
        # energy of defect complex
        binding_energy = stable_charges[name][1]
        # subtracting sum of energies of single defects
        binding_energy += -1 * sum([abs(v) * stable_charges[single_defects[el]][1]
                                    for el,v in defect_complex.delta_atoms.items()])
    
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
        npoints = 1000
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
    
    
    def defect_concentrations(self, chemical_potentials, temperature=300, fermi_level=0.):
        """
        Give list of all concentrations at specified efermi
        args:
            chemical_potentials = {Element: number} is a dictionary of chemical potentials  
            temperature = temperature to produce concentrations from
            fermi_level: (float) is fermi level relative to valence band maximum
                Default efermi = 0 = VBM energy
        returns:
            list of dictionaries of defect concentrations
        """
        concentrations = []
        for dfct in self.entries:
            concentrations.append({
                'conc':
                dfct.defect_concentration(
                    self.vbm, chemical_potentials=chemical_potentials, temperature=temperature, fermi_level=fermi_level),
                'name':
                dfct.name,
                'charge':
                dfct.charge
            })

        return concentrations
    

    def defect_concentrations_stable_charges(self, chemical_potentials, temperature=300, fermi_level=0.):
        """
        Give list of concentrations of defects in their stable charge state at a specified efermi.

        Parameters
        ----------
        chemical_potentials : (Dict)
            Dictionary of chemical potentials.
        temperature : (float), optional
            Temperature. The default is 300.
        fermi_level : (float), optional
            Position of the Fermi level. The default is 0 (vbm).

        Returns
        -------
        conc_stable : (list)
            List of dictionaries with concentrations of defects with stable charge states at a given efermi.
        """
        concentrations = self.defect_concentrations(chemical_potentials,temperature=temperature,fermi_level=fermi_level)
        stable_charges = self.stable_charges(chemical_potentials,fermi_level=fermi_level)
        conc_stable = []
        for c in concentrations:
            n = c['name']
            q = c['charge']
            for name in stable_charges:
                charge = stable_charges[name][0]
                if n == name and q == charge:
                    conc_stable.append(c)
        return conc_stable


    def defect_concentrations_total(self, chemical_potentials, temperature=300, fermi_level=0.):
        """
        Calculate the sum of the defect concentrations in every charge state for every defect specie.

        Parameters
        ----------
        chemical_potentials : (Dict)
            Dictionary of chemical potentials.
        temperature : (float), optional
            Temperature. The default is 300.
        fermi_level : (float), optional
            Position of the Fermi level. The default is 0 (vbm).

        Returns
        -------
        total_concentrations : (Dict)
            Dictionary with names of the defect species as keys and total concentrations as values.

        """
        
        total_concentrations = {}
        for name in self.names:
            total_concentrations[name] = 0
            for d in self.defect_concentrations(chemical_potentials=chemical_potentials,
                                                temperature=temperature,fermi_level=fermi_level):
                if d['name'] == name:
                    total_concentrations[name] += d['conc']
        
        return total_concentrations
                    
            
    def equilibrium_fermi_level(self, chemical_potentials, bulk_dos, temperature = 300):
        """
        Solve for the Fermi energy self-consistently as a function of T
        Observations are Defect concentrations, electron and hole conc
        Args:
            temperature: Temperature to equilibrate fermi energies for
            chemical_potentials: dict of chemical potentials to use for calculation fermi level
            bulk_dos: bulk system dos (pymatgen Dos object)
        Returns:
            Fermi energy dictated by charge neutrality
        """

        fdos = FermiDos(bulk_dos, bandgap=self.band_gap)
        _,fdos_vbm = fdos.get_cbm_vbm()

        def _get_total_q(ef):

            qd_tot = sum([
                d['charge'] * d['conc']
                for d in self.defect_concentrations(
                    chemical_potentials=chemical_potentials, temperature=temperature, fermi_level=ef)
            ])
            qd_tot += fdos.get_doping(fermi_level=ef + fdos_vbm, temperature=temperature)
            return qd_tot

        return bisect(_get_total_q, -1., self.band_gap + 1.)
            

    def _get_frozen_charge(self,frozen_defect_concentrations, chemical_potentials, temperature=300, fermi_level=0.):
        
        conc = self.defect_concentrations(chemical_potentials,temperature,fermi_level)
        total_conc = self.defect_concentrations_total(chemical_potentials,temperature,fermi_level)
        frozen_conc = {fr['name']:fr['conc'] for fr in frozen_defect_concentrations}
        
        qd_tot = 0
        for d in conc:
            name = d['name'] 
            entry = self.find_entry_by_name_and_charge(name,d['charge'])
            c = d['conc']
            #handle defect complex case
            if entry.__class__.__name__ == 'DefectComplexEntry': 
                complex_elements = [d.site.specie.symbol for d in entry.defect_list]
                for el in complex_elements:
                    for name in frozen_conc.keys():
                        if el in name:
                            if name in total_conc.keys():
                                c = c * (frozen_conc[name] / total_conc[name])
                            else:
                                print(f'Warning: Defect with name {name} is not in defect entries and will not affect the calculation')
                qd_tot += d['charge'] * c
            
            #handle single defect case
            else: 
                # D1 group
                if name in frozen_conc.keys():
                    if name in total_conc.keys():
                        qd_tot += d['charge'] * d['conc'] * (frozen_conc[name] / total_conc[name])
                    else:
                        print(f'Warning: Defect with name {name} is not in defect entries and will not affect the calculation')
                # D2 group
                else: 
                    qd_tot += d['charge'] * d['conc']
                
            return qd_tot
            
        
    def non_equilibrium_fermi_level(self, frozen_defect_concentrations, chemical_potentials, bulk_dos, 
                                        external_defects=[], temperature=300):
        """
        Solve charge neutrality in non-equilibrium conditions. The contribution to the total charge concentration
        of the defects can arise from 3 different contributions (groups):
            - D1 : frozen defects with defect specie present in defect entries
            - D2 : normal defects with defect specie present in defect entries
            - D3 : external defects not present in defect entries, with fixed concentration and charge
            
        The group D1 is the less straightforward. Frozen defects can also be relative to another system and another temperature,
        but only the ones with names found in the defect entries will be computed. 
        This part will account for a fixed installed distribution of intrisic defects that are allowed to 
        change their charge state. One example would be to equilibrate the system with concentrations
        of intrinsic defects installed in another phase at a higher temperature, that are considered constant 
        (kinetically trapped) as the temperature lowers and the phases changes. The variation of the charge states is 
        accounted for in the solution of the charge neutrality by defining a charge state distribution for every defect
        specie, which depends of the fermi level. 
        If a defect entry is not found in group D1 is considered to belong to group D2. The number of elements of 
        D1 U D2 will be equal to the muber of defect entries.
        The charge concentration associated to group D2 is treated as in the "equilibrium_fermi_level" function.
        The charge concentration associated to group D3 is a number and thus not Fermi level dependent.
        
        Parameters
        ----------
        frozen_defect_concentrations : (list)
            List of defect concentrations. Most likely generated with the defect_concentrations() method. It is not
            recommended to generate this manually.
        chemical_potentials : (Dict)
            Dictionary of chemical potentials in the format {Element('el'):chempot}.
        bulk_dos : (CompleteDos object)
            Pymatgen CompleteDos object of the DOS of the bulk system.
        external_defects : (list)
            List of external defect concentrations (not present in defect entries).
        temperature : (float), optional
            Temperature to equilibrate the system to. The default is 300.

        Returns
        -------
        (float)
            Fermi level dictated by charge neutrality .
        """
        
        fdos = FermiDosCarriersInfo(bulk_dos, bandgap=self.band_gap)
        _,fdos_vbm = fdos.get_cbm_vbm()

        def _get_total_q(ef):
            
            # get groups D1 and D2
            qd_tot = self._get_frozen_charge(frozen_defect_concentrations, chemical_potentials,temperature,fermi_level=ef)

            #external fixed defects - D3
            for d_ext in external_defects:
                qd_tot += d_ext['charge']*d_ext['conc']
                
            qd_tot += fdos.get_doping(fermi_level=ef + fdos_vbm, temperature=temperature)

            return qd_tot
                       
        return bisect(_get_total_q, -1., self.band_gap + 1.)

    
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
            
    
    def plot(self,mu_elts=None,xlim=None, ylim=None, title=None, fermi_level=None, 
             plotsize=1, fontsize=1.2, show_legend=True, format_legend=False, 
             order_legend=False, get_subplot=False, subplot_settings=None):
        """
        Produce defect Formation energy vs Fermi energy plot
        Args:
            mu_elts:
                a dictionnary of {Element:value} giving the chemical
                potential of each element
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
        
        # creating a dictionary with names of single defect entry as keys and
        # a list of tuples (charge,formation_energy) as values
        # Every defect name identifies a type of defect
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
     

    def _get_formatted_legend(self,name):
        
        if '-' not in [c for c in name]:        
            flds = name.split('_')
            if 'Vac' == flds[0]:
                base = '$V'
                sub_str = '_{' + flds[1] + '}$'
            elif 'Sub' == flds[0]:
                flds = name.split('_')
                base = '$' + flds[1]
                sub_str = '_{' + flds[3] + '}$'
            elif 'Int' == flds[0]:
                base = '$' + flds[1]
                sub_str = '_{int}$'
            else:
                base = name
                sub_str = ''
    
            return  base + sub_str
        
        else:
            label = ''
            names = name.split('-')
            for name in names:
                flds = name.split('_')
                if '-' not in flds:
                    if 'Vac' == flds[0]:
                        base = '$V'
                        sub_str = '_{' + flds[1] + '}$'
                    elif 'Sub' == flds[0]:
                        flds = name.split('_')
                        base = '$' + flds[1]
                        sub_str = '_{' + flds[3] + '}$'
                    elif 'Int' == flds[0]:
                        base = '$' + flds[1]
                        sub_str = '_{int}$'
                    else:
                        base = name
                        sub_str = ''
            
                    if names.index(name) != (len(names) - 1):
                        label += base + sub_str + '-'
                    else:
                        label += base + sub_str
            
            return label
    
    
    def plot_binding_energies(self, names, xlim=None, ylim=None, size=1):
        """
        Plot binding energies for complex of defects as a function of the fermi level
        Args:
            names: 
                List of strings with names of defect complex as assigned in SingleDefectData object
            xlim:
                Tuple (min,max) giving the range of the x (fermi energy) axis
            ylim:
                Tuple (min,max) giving the range for the formation energy axis
            size:
                Float multiplier to change plot size
        """
        
        plt.figure(figsize=(8*size,6*size))
        matplotlib.rcParams.update({'font.size': 10*1.8*size}) 
        
        # BINDING ENERGY
        if xlim==None:
            xlim = (-0.5,self.band_gap+0.5)
        # building array for x values (fermi level)    
        ef = np.arange(xlim[0],xlim[1]+0.1,(xlim[1]-xlim[0])/200)        
        binding_energy = np.zeros(len(ef))
        
        # getting binding energy at different fermi levels for every name in list
        for name in names:
            for i in range(0,len(ef)):
                binding_energy[i] = self.binding_energy(name,fermi_level=ef[i])            
            plt.plot(ef,binding_energy, linewidth=2.5*size,label=name)
            
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
                
        """
        
        plt.figure(figsize=(10*size,10*size)) 
        
        if ylim == None:
            ylim = (-0.5,self.band_gap +0.5)
        
        # get all charge transition levels
        charge_transition_levels = self.charge_transition_levels()
        
        # number of defect names
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
                flds = name.split('_')
                if 'Vac' == flds[0]:
                    base = '$V'
                    sub_str = '_{' + flds[1] + '}$'
                elif 'Sub' == flds[0]:
                    flds = name.split('_')
                    base = '$' + flds[1]
                    sub_str = '_{' + flds[3] + '}$'
                elif 'Int' == flds[0]:
                    base = '$' + flds[1]
                    sub_str = '_{int}$'
                else:
                    base = name
                    sub_str = ''
            
                x_ticks_labels[x_ticks_labels.index(name)] = base + sub_str
        
        
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
        
            
            
            
        
        
        
        
        
        
        
        
        
        
        