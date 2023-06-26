
import warnings
import json
import os.path as op
import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt
from monty.json import MSONable, MontyEncoder
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram, PDPlotter
from pymatgen.symmetry.groups import SpaceGroup
from pynter.tools.format import format_composition
import copy


class Chempots(MSONable):
    
    def __init__(self,chempots_dict,ndecimals=6):
        """
        Class to handle set of chemical potentials. Behaves like a python dictionary.
        The dictionary needs to be set with element symbols as keys and chemical potentials 
        as values.

        Parameters
        ----------
        chempots_dict : (dict)
            Dictionary of chemical potentials in the format {el:value}.
        ndecimals : (int)
            Round the chemical potentials to this number of decimals. If None
            the numbers are left untouched. the default is 6.
            
        """
        if ndecimals:
            self._mu = {el:round(v,ndecimals) for el,v in chempots_dict.items()}
        else:
            self._mu = chempots_dict
    
    @property
    def mu(self):
        """
        Dictionary with chemical potentials
        """
        return self._mu

    def __str__(self):
        return self.mu.__str__()
    
    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self.mu)

    def __iter__(self):
        return self.mu.keys().__iter__()

    def __getitem__(self,el):
        return self.mu[el]

    def __setitem__(self,el,value):
        self.mu[el] = value
        return
    
    def __delitem__(self,el):
        del self.mu[el]
        return
    
    def __eq__(self, other):
        if isinstance(other, dict):
            return self.mu == other
        elif isinstance(other, Chempots):
            return self.mu == other.mu
        else:
            return False

    def keys(self):
        return self.mu.keys()

    def values(self):
        return self.mu.values()
    
    def items(self):
        return self.mu.items()

    def copy(self):
        return Chempots(copy.deepcopy(self.mu))
    
    def update(self, other):
        if isinstance(other, dict):
            for key, value in other.items():
                self.mu[key] = value
        else:
            for key, value in other:
                self.mu[key] = value
    
    
    def as_dict(self):
        """
        Json-serializable dict representation of a Chempots object. 
        
        Returns
        -------
        dict
            Json-serializable dict of a Chempots object.
        """
        d = {
        "@module": self.__class__.__module__,
        "@class": self.__class__.__name__,
        "chempots":self.mu
        }
        return d
    
    @staticmethod
    def from_dict(d):
        """
        Constructor of Chempots object from dictionary representation.
        
        Parameters
        ----------
        d : dict
        
        Returns
        -------
        Chempots object.
        """
        chempots = d["chempots"]
        return Chempots(chempots)
    
    @staticmethod
    def from_pmg_elements(d):
        """
        Construct a Chempots object from a dictionary in the format {Element:value}. Useful to
        interface with pymatgen which uses Element class to address elements instead of just their symbols.

        Parameters
        ----------
        d : (dict)
            Dictionary in the format {Element:value}.

        Returns
        -------
        Chempots.
        """
        return Chempots({el.symbol:value for el,value in d.items()})
    
    def to_pmg_elements(self):
        """
        Convert the Chempots object to a dictionary in the format {Element:value}. 
        Useful to interface with pymatgen which uses Element class to address elements 
        instead of just their symbols.

        Returns
        -------
        Dictionary in the format {Element:value}
        """
        return {Element(el):value for el,value in self.mu.items()}
    
    
    def get_absolute(self,mu_refs):
        """
        Get Chempots object with chemical potentials converted to absolute values (mu + mu_ref)

        Parameters
        ----------
        mu_refs : (dict or Chempots)
            Chempots object or just dictionary with elemental chemical potentials (mu0).

        Returns
        -------
        Chempots values converted to absolute (mu + mu_ref).
        """
        return Chempots({el:self.mu[el] + mu_refs[el] for el in self.mu})
        
    
    def get_referenced(self,mu_refs):
        """
        Get Chempots object with chemical potentials converted to referenced values (mu - mu_ref)

        Parameters
        ----------
        mu_refs : (dict or Chempots)
            Chempots object or just dictionary with elemental chemical potentials (mu0).

        Returns
        -------
        Chempots values converted to referenced (mu - mu_ref).
        """
        return Chempots({el:self.mu[el] - mu_refs[el] for el in self.mu})
            
        
        
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
        self.pd = phase_diagram if phase_diagram else None
        if mu_refs:
            self.mu_refs = mu_refs
        elif self.pd:
            self.mu_refs = PDHandler(self.pd).get_chempots_reference()
        else:
            self.mu_refs = mu_refs
            warnings.warn('Neither PhaseDiagram or reference chempots have been provided, conversions btw ref and abs value will not be possible',UserWarning)
        self.are_chempots_delta = are_chempots_delta


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
        self.are_chempots_delta = False
        return

        
    def set_to_referenced(self):
        """
        Set reservoir dictionary to referenced values of chempots.
        """
        new_res_dict = self.get_referenced_res_dict()
        self.res_dict = new_res_dict
        self.are_chempots_delta = True
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
                    label = '$\Delta \mu_{\text{%s}}$' %el
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
    
    
    
    
class PDHandler:
    
    def __init__(self,phase_diagram):
        """
        Class to generate and handle Pymatgen phase diagram more rapidly.

        Parameters
        ----------
        phase_diagram : (PhaseDiagram)
            Pymatgen PhaseDiagram object
        """
        self.pd = phase_diagram 
        self.mu_refs = self.get_chempots_reference()


    def calculate_single_chempot(self,comp,chempots_ref):
        """
        Calculate referenced chemical potential in a given composition and given
        the chemical potentials of the other elements.

        Parameters
        ----------
        comp : (Pymatgen Composition object)
            Compositions of the phase.
        chempots_ref : (Dict)
            Dictionary with element symbols as keys and respective chemical potential as value ({el:chempot}).
            The chemical potentials used here are the ones relative to the reference (delta_mu)

        Returns
        -------
        mu : (float)
            Chemical potential.
        """        
        form_energy = self.get_formation_energy_from_stable_comp(comp)
        mu = form_energy
        for el,coeff in comp.items():
            if el.symbol in chempots_ref:
                mu += -1*coeff * chempots_ref[el.symbol]
            else:
                factor = coeff
        
        return mu/factor


    def get_all_boundaries_chempots(self,comp):
        """
        Get chemical potentials at the corners of the stability ta given composition.

        Parameters
        ----------
        comp : 
            Pymatgen Composition object.

        Returns
        -------
        chempots : (Chempots)
            Chempots object.
        """
        chempots_pmg = self.pd.get_all_chempots(comp)
        for r,mu in chempots_pmg.items():
            chempots_pmg[r] = {k:v.item() for k,v in mu.items()} # convert from np.float64
        chempots = {r:Chempots.from_pmg_elements(mu) for r,mu in chempots_pmg.items()}
        return chempots

        
    def get_chempots_reference(self):
        """
        Gets elemental reference compounds and respective e.p.a with Pymatgen el_ref attribute in PhaseDiagram class.

        Returns
        -------
        Chempots object
        """
        chempots_ref = {}
        pd = self.pd
        for el in pd.el_refs:
            chempots_ref[el] = pd.el_refs[el].energy_per_atom
        chempots_ref = {k: v for k, v in sorted(chempots_ref.items(), key=lambda item: item[0])}
        return Chempots.from_pmg_elements(chempots_ref)
       

    def get_dataframe(self):
        """
        Generate pandas DataFrame with columns 'Composition, Structure, Formation energy'.
        To display a string for 'Structure' the entry needs to be a ComputedStructureEntry (see pymatgen docs).
        """
        phases = []
        for e in self.pd.stable_entries:
            d = {}
            d['Composition'] = format_composition(e.composition.reduced_formula)
            sg = SpaceGroup(e.structure.get_space_group_info()[0])
            crystal_system, sg_string = sg.crystal_system , sg.to_latex_string()
            if e.__class__.__name__ == 'ComputedStructureEntry':
                d['Structure'] = f'{crystal_system.capitalize()} ({sg_string})'
            else:
                d['Structure'] = None
            d['Formation energy p.a (eV)'] = np.around(self.pd.get_form_energy_per_atom(e),decimals=2)
            phases.append(d)
        df = DataFrame(phases)
        return df

        
    def get_entries_from_comp(self,comp):
       """
       Get a list of entries corrisponding to the target composition.
       
       Parameters
       ----------
       comp : (Pymatgen Composition object)
       Returns
       -------
       List of Pymatgen PDEntry objects
       """
       target_entries=[]
       pd = self.pd
       for e in pd.all_entries:
           if e.composition.reduced_composition == comp:
               target_entries.append(e)
       if target_entries != []:
           return target_entries
       else:
           raise ValueError('No entry has been found for target composition:%s' %comp.reduced_formula)
               

    def get_formation_energies_from_comp(self,comp):
        """
        Get dictionary of formation energies for all entries of a given reduced composition

        Parameters
        ----------
        comp : (Pymatgen Composition object)
        Returns
        -------
        from_energies : (Dict)
            Dictionary with PDEntry objects as keys and formation energies as values.
        """
        pd = self.pd
        form_energies = {}
        for e in self.get_entries_from_comp(comp):
            factor = e.composition.get_reduced_composition_and_factor()[1]
            form_energies[e] = pd.get_form_energy(e)/factor
        return form_energies
            

    def get_formation_energy_from_stable_comp(self,comp):
        """
        Get formation energy of a target stable composition

        Parameters
        ----------
        comp : (Pymatgen Composition object)

        Returns
        -------
        Formation energy (float)
        """
        pd = self.pd
        entry = self.get_stable_entry_from_comp(comp)
        factor = entry.composition.get_reduced_composition_and_factor()[1]
        return pd.get_form_energy(entry)/factor
        

    def get_phase_boundaries_chempots(self,comp,chempot_ref):
        """
        Given a composition and a fixed chemical potential, this function analises the composition of the boundary phases
        and the associated chemical potentials at the boundaries. Only works for 3 component PD.

        Parameters
        ----------
        comp : (Pymatgen Composition object)
            Composition of the phase you want to get the chemical potentials at the boundary.
        chempot_ref : (Dict)
            Dictionary with fixed element symbol as key and respective chemical potential as value ({el:chempot}).
            The chemical potential here is the referenced value. 
        Returns
        -------
        chempots : (Dict)
            Dictionary with compositions at the boundaries as keys and delta chemical potentials as value.
        """       
        chempots = {}
        comp1,comp2 = self.get_phase_boundaries_compositions(comp, chempot_ref)
        
        boundary = '-'.join([comp1.reduced_formula,comp.reduced_formula])
        chempots[boundary] = self.solve_phase_boundary_chempots(comp1, comp, chempot_ref)
        
        boundary = '-'.join([comp.reduced_formula,comp2.reduced_formula])
        chempots[boundary] = self.solve_phase_boundary_chempots(comp, comp2, chempot_ref)
        
        return chempots      
                  
    
    def get_phase_boundaries_compositions(self,comp,chempot_ref):
        """
        Get compositions of phases in boundary of stability with a target composition given a fixed chemical potential 
        on one component. Currently only works for 3-component PD (to check). 
        Used Pymatgen GrandPotentialPhaseDiagram class. The fixed chemical potential is the referenced value that is
        converted in the global value for the analysis with the GrandPotentialPhaseDiagram class.

        Parameters
        ----------
        comp : (Pymatgen Composition object)
            Target composition for which you want to get the bounday phases.
        chempot_ref : (Dict)
            Dictionary with fixed element symbol as key and respective chemical potential as value ({el:chempot}).
            The chemical potential is the referenced value

        Returns
        -------
        comp1,comp2 : (Pymatgen Composition object)
            Compositions of the boundary phases given a fixed chemical potential for one element.
        """
        chempot_ref = Chempots(chempot_ref)
        fixed_chempot = chempot_ref.get_absolute(self.mu_refs).to_pmg_elements()
        
        entries = self.pd.all_entries
        gpd = GrandPotentialPhaseDiagram(entries, fixed_chempot)
        stable_entries = gpd.stable_entries
        comp_in_stable_entries = False
        for e in stable_entries:
            if e.original_comp.reduced_composition == comp:
                comp_in_stable_entries = True
        if comp_in_stable_entries == False:
            raise ValueError('Target composition %s is not a stable entry for fixed chemical potential: %s' %(comp.reduced_formula,fixed_chempot))
        
        el = comp.elements[0]
        x_target = comp.get_wt_fraction(el)
        x_max_left = 0
        x_min_right = 1
        for e in stable_entries:
            c = e.original_comp
            if c != comp:
                x = c.get_wt_fraction(el)
                if x < x_target and x >= x_max_left:
                    x_max_left = x
                    comp1 = c
                if x > x_target and x <= x_min_right:
                    x_min_right = x
                    comp2 = c
        return comp1.reduced_composition,comp2.reduced_composition   
                  
                    
    def solve_phase_boundary_chempots(self,comp1,comp2,chempot_ref):
        """
        Given a fixed chemical potential, gets the values of the remaining two chemical potentials
        in the boundary between two phases (region where the two phases coexist). Only works for 3-component PD (to check).
        Given a phase P1 (formula AxByOz) and a phase P2 (formula AiBjOk) the chemical potentials have to satisfy the conditions:
            form_energy(P1) = x*mu(A) + y*mu(B) +z*mu(O)
            form_energy(P2) = i*mu(A) + j*mu(B) +k*mu(O)
        From these conditions the values of mu(A) and mu(B) are determined given a fixed value of mu(O).
        All of the chemical potentials used here are delta_mu, i.e. relative to the elemental phase(delta_mu(O) = mu(O) - mu_ref(O))

        Parameters
        ----------
        comp1,comp2 : (Pymatgen Composition object)
            Compositions of the two phases at the boundary.
        chempot_ref : (Dict)
            Dictionary with fixed element symbol as key and respective chemical potential as value ({el:chempot}). The chemical potential
            used here is the one relative to the reference (delta_mu)

        Returns
        -------
        chempots_boundary : (Dict)
            Dictionary of chemical potentials.
        """
        chempots_boundary ={}
        chempot_ref = Chempots(chempot_ref)
        mu_fixed = chempot_ref.to_pmg_elements()
        for el,chempot in mu_fixed.items():
            el_fixed, mu_fixed = el, chempot
        e1 = self.get_formation_energy_from_stable_comp(comp1)
        e2 = self.get_formation_energy_from_stable_comp(comp2)
        
        coeff1 = []
        coeff2 = []
        # order of variables (mu) will follow the order of self.chempots_reference which is alphabetically ordered
        mu_refs = self.mu_refs.to_pmg_elements()
        for el in mu_refs:
            if el != el_fixed:
                if el not in comp1:
                    coeff1.append(0)
                else:
                    coeff1.append(comp1[el])
                if el not in comp2:
                    coeff2.append(0)
                else:
                    coeff2.append(comp2[el])                
        a = np.array([coeff1,coeff2])
        
        const1 = e1 - comp1[el_fixed]*mu_fixed if el_fixed in comp1 else e1 
        const2 = e2 - comp2[el_fixed]*mu_fixed if el_fixed in comp2 else e2
        b = np.array([const1,const2])
        x = np.linalg.solve(a, b)
        # output will follow the order given in input
        counter = 0
        for el in mu_refs:
            if el != el_fixed:
                chempots_boundary[el] = x[counter]
                counter += 1
        chempots_boundary[el_fixed] = mu_fixed        
        
        return Chempots.from_pmg_elements(chempots_boundary) 


    def get_plot(self,**kwargs):
        """
        Get plot with Pymatgen
        """
        PDPlotter(self.pd,show_unstable=0,backend='matplotlib').get_plot(**kwargs)
        return plt


    def get_stability_diagram(self,elements,size=None):
        """
        Method to get stability diagram with 'get_chempot_range_map_plot' method in pymatgen

        Parameters
        ----------
        elements : (list)
            List with strings of the elements to be used as free variables.
        size : (tuple)
            New size in inches.

        Returns
        -------
        plt : 
            Matplotlib object
        """
        pd = self.pd
        elements = [Element(el) for el in elements]
        PDPlotter(pd).get_chempot_range_map_plot(elements)
        if size:
            fig = plt.gcf()
            fig.set_size_inches(size[0],size[1])
        
        return plt


    def get_stable_entry_from_comp(self,comp):
        """
        Get the PDEntry of the stable entry of a target composition

        Parameters
        ----------
        comp : (Pymatgen Composition object)
        Returns
        -------
        Pymatgen PDEntry object
        """
        target_entry=None
        pd = self.pd
        for e in pd.stable_entries:
            if e.composition.reduced_composition == comp:
                target_entry = e
                break
        if target_entry is not None:
            return target_entry
        else:
            raise ValueError('No stable entry has been found for target composition:%s' %comp.reduced_formula)

    
