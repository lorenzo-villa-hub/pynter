
import numpy as np
from monty.json import MSONable
from pymatgen.core.periodic_table import Element
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



            
    

    
