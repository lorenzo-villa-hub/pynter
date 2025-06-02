
import numpy as np
from itertools import product
from monty.json import MSONable
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
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


def chempot_ideal_gas(mu0, temperature,partial_pressure):
    """
    Get chemical potential at a given temperature and partial pressure. The chemical potential in standard conditions (mu0)
    has to be know.

    Parameters
    ----------
    mu0 : (float)
        Value of chemical potential in standard conditions.
    temperature : (float)
        Temperature. in Kelvin
    partial_pressure : (float)
        Partial pressure.

    Returns
    -------
    chempot : (float)
        Value of chemical potential at given T and p/p0.
    """
    kb = 8.6173324e-5  # eV / K
    chempot = mu0 + 0.5*kb*temperature*np.log(partial_pressure)
    return chempot


def barycenter_chemical_potentials_absolute(composition,
                                            energy,
                                            oxygen_chempot_absolute,
                                            mu_refs,
                                            min_absolute_chempots=None,
                                            max_absolute_chempots=None):
    """
    Compute the barycenter of the feasible region for relative chemical potentials,
    constrained by:
      - fixed absolute chemical potential of oxygen
      - Total energy of target phase
      - lower and/or upper chemical potential limits for each element
    
    Parameters:
    -----------
    composition: str or pymatgen.core.Composition
        Target composition.
    formation_energy: float
        Formation energy (eV/f.u.)
    oxygen_mu_relative: float
        fixed chemical potential of oxygen relative to the oxygen molecule.
    min_relative_chempots: dict or Chempots
        Lower limit of chemical potentials, relative values ({element:chempot}) 
    max_relative_chempots: dict or Chempots
        Higher limit of chemical potentials, relative values ({element:chempot})
    
    Returns:
        Dictionary with chemical potentials, taken from the center of the
        allowed N-1 dimensional hyperplane.
    """ 
    if isinstance(composition, str):
        composition = Composition(composition)

    formation_energy = energy - sum([number*mu_refs[el.symbol] for el,number in composition.items()])
    oxygen_chempot_relative = oxygen_chempot_absolute - mu_refs['O']
    min_relative_chempots = Chempots(min_absolute_chempots).get_referenced(mu_refs) if min_absolute_chempots else None
    max_absolute_chempots = Chempots(max_absolute_chempots).get_referenced(mu_refs) if max_absolute_chempots else None
    
    barycenter_chempots_relative = barycenter_chemical_potentials_relative(composition=composition,
                                                   formation_energy=formation_energy,
                                                   oxygen_chempot_relative=oxygen_chempot_relative,
                                                   min_relative_chempots=min_relative_chempots,
                                                   max_relative_chempots=max_absolute_chempots)
    
    return Chempots(barycenter_chempots_relative).get_absolute(mu_refs)
        

    

def barycenter_chemical_potentials_relative(composition,
                                           formation_energy,
                                           oxygen_chempot_relative,
                                           min_relative_chempots=None,
                                           max_relative_chempots=None):
    """
    Compute the barycenter of the feasible region for relative chemical potentials,
    constrained by:
      - fixed relative chemical potential of oxygen
      - Formation energy of target phase
      - lower and/or upper chemical potential limits for each element
    
    Parameters:
    -----------
    composition: str or pymatgen.core.Composition
        Target composition.
    formation_energy: float
        Formation energy (eV/f.u.)
    oxygen_mu_relative: float
        fixed chemical potential of oxygen relative to the oxygen molecule.
    min_relative_chempots: dict or Chempots
        Lower limit of chemical potentials, relative values ({element:chempot}) 
    max_relative_chempots: dict or Chempots
        Higher limit of chemical potentials, relative values ({element:chempot})
    
    Returns:
        Dictionary with chemical potentials, taken from the center of the
        allowed N-1 dimensional hyperplane.
    """
    muO_relative = oxygen_chempot_relative
    if isinstance(composition, str):
        composition = Composition(composition)
    
    el_amount_dict = composition.get_el_amt_dict()
    if "O" not in el_amount_dict:
        raise ValueError("Oxygen must be present in the composition.")

    n_O = el_amount_dict['O']
    non_oxygen_elements = [el for el in el_amount_dict if el !='O']
    constraint = formation_energy - n_O * muO_relative
    chempot_boundaries = {}
    for el in non_oxygen_elements:
        if min_relative_chempots and el in min_relative_chempots.keys():
            mu_min = min_relative_chempots[el]
        else:
            mu_min = float("-inf")
        if max_relative_chempots and el in max_relative_chempots.keys():
            mu_max = max_relative_chempots[el]
        else:
            mu_max = 0.0
        chempot_boundaries[el] = (mu_min,mu_max)
    
    allowed_vertices = []
    # iterate over all possible permutations of chempot ranges
    for node in product(*[chempot_boundaries[el] for el in non_oxygen_elements]):
        for idx in range(len(non_oxygen_elements)):
            mu_values = list(node)
            el_free = non_oxygen_elements[idx]
            n_el_free = el_amount_dict[el_free]
            
            fixed_sum = sum(el_amount_dict[non_oxygen_elements[i]] * mu_values[i]
                            for i in range(len(non_oxygen_elements)) if i != idx)
            
            mu_free = (constraint - fixed_sum) / n_el_free  # solve free chemical potentials
            
            mu_min, mu_max = chempot_boundaries[el_free]
            if mu_min <= mu_free <= mu_max:  #if solved chempot is within limits, we add vertex to N-1 chempot hyperplane
                mu_values[idx] = mu_free
                vertex = dict(zip(non_oxygen_elements, mu_values))
                allowed_vertices.append(vertex)
                
    if not allowed_vertices:
        raise ValueError("No feasible chemical potential points found under the given constraints.")
        
    barycenter = {el:0 for el in non_oxygen_elements}
    for vertex in allowed_vertices:
        for el, mu in vertex.items():
            barycenter[el] += mu
    for el in barycenter.keys():
        barycenter[el] /= len(allowed_vertices)
    barycenter['O'] = muO_relative
    
    return barycenter
            
    

    
