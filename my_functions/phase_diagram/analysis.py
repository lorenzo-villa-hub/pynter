
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram


class ChempotAnalysis:
    
    def __init__(self,computed_phases):
        """
        Initializes class to analyse chemical potentials of a phase diagram generated with Pymatgen
        
        Parameters
        ----------
        computed_phases : (Dict)
            Dictionary with strings with phase reduced formula as keys and energy per formula unit as values.
            The strings with reduced formula are trasformed in Pymatgen Composition objects to be used in this class.
        """
        self._computed_phases = {Composition(phase):computed_phases[phase] for phase in computed_phases}
        self._chempots_reference = self._get_chempots_reference()


    @property
    def computed_phases(self):
        return self._computed_phases


    @property
    def chempots_reference(self):
        return self._chempots_reference
    
    
    # def _get_chempots_reference(self):
    #     # gets reference chempots from e.p.a of elemental compounds in self.computed_phases
    #     chempots_ref = {}
    #     for comp in self.computed_phases:
    #             if comp.is_element:
    #                 chempots_ref[comp.elements[0]] = self.computed_phases[comp]/comp.num_atoms            
    #     return chempots_ref
    

    def _get_chempots_reference(self):
        # gets elemental reference compounds and respective e.p.a with Pymatgen el_ref attribute in PhaseDiagram class
        chempots_ref = {}
        pd = PhaseDiagram(self.get_pd_entries())
        for el in pd.el_refs:
            chempots_ref[el] = pd.el_refs[el].energy_per_atom
        return chempots_ref
   
    def get_pd_entries(self):
        """
        Build list of PDEntry object used by Pymatgen to generate PhaseDiagram.
        Returns
        -------
        List of PDEntries
        """  
        entries = []
        for comp in self.computed_phases:
            entry = PDEntry(comp,self.computed_phases[comp])
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
                    
            
  #  def get_chempots_limit(self,comp1,comp2):
        
                  
                    
        