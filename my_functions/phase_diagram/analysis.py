
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PDEntry

class ChempotAnalysis:
    
    def __init__(self,computed_phases):       
        self._computed_phases = computed_phases
        self._chempots_reference = self._get_chempots_reference()


    @property
    def computed_phases(self):
        return self._computed_phases


    @property
    def chempots_reference(self):
        return self._chempots_reference
    
    
    def _get_chempots_reference(self):
        # gets reference chempots from e.p.a of elemental compounds in self.computed_phases
        chempots_ref = {}
        for phase in self.computed_phases:
                if len(Composition(phase).elements)==1:
                    el_ref = phase
                    comp_ref = Composition(el_ref)
                    chempots_ref[comp_ref.elements[0]] = self.computed_phases[el_ref]/comp_ref.num_atoms       
        return chempots_ref
    
    
    def get_pd_entries(self):
        """
        Build list of PDEntry object used by Pymatgen to generate PhaseDiagram
        Returns
        -------
        List of PDEntries
        """  
        entries = []
        for phase in self.computed_phases:
            comp = Composition(phase)
            entry = PDEntry(comp,self.computed_phases[phase])
            entries.append(entry)           
        return entries
 
    
    def get_chempots_delta(self,chempots_abs):
        """
        Subtract energy per atom of elemental phases to dictionary of chemical potentials ({Element:chempot})
        The energy of the elemental phase is taken from self.computed_phases
        Parameters
        ----------
        chempots_abs : (Dict)
            Dictionary of absolute chemical potentials, format ({Pymatgen Element object:chempot value}).
        Returns
        -------
        chempots_delta : (Dict)
            Dictionary of chemical potentials relative to elemental phase, format ({Pymatgen Element object:chempot value}).
        """
        import copy
        chempots_delta = copy.deepcopy(chempots_abs)
        for el in chempots_abs:
            chempots_delta[el] += (-1)*self.chempots_reference[el]            
        return chempots_delta
                    
            
                    
                    
        