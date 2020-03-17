
import numpy as np
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
    

    def _get_chempots_reference(self):
        # gets elemental reference compounds and respective e.p.a with Pymatgen el_ref attribute in PhaseDiagram class
        # the difference with my version is that this considers having more than 1 calculation and takes the minimum e.p.a
        chempots_ref = {}
        pd = PhaseDiagram(self.get_pd_entries())
        for el in pd.el_refs:
            chempots_ref[el] = pd.el_refs[el].energy_per_atom
        chempots_ref = {k: v for k, v in sorted(chempots_ref.items(), key=lambda item: item[0])}
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
                    
            
    def get_chempots_boundary(self,comp1,comp2,fixed_chempot):
        """
        Given a fixed chemical potential, gets the values of the remaining two chemical potentials
        in the boundary between two phases (region where the two phases coexist). Only works for 3-component PD (to check).
        Given a phase P1 (formula AxByOz) and a phase P2 (formula AiBjOk) the chemical potentials have to satisfy the conditions:
            form_energy(P1) = x*mu(A) + y*mu(B) +z*mu(O)
            from_energy(P2) = i*mu(A) + j*mu(B) +k*mu(O)
        From these conditions the values of mu(A) and mu(B) are determined given a fixed value of mu(O).
        All of the chemical potentials used here are delta_mu, i.e. relative to the elemental phase(delta_mu(O) = mu(O) - mu_ref(O))

        Parameters
        ----------
        comp1,comp2 : (Pymatgen Composition object)
            Compositions of the two phases at the boundary.
        fixed_chempot : (Dict)
            Dictionary with fixed Element as key and respective chemical potential as value ({Element:chempot}).

        Returns
        -------
        chempots_boundary : (Dict)
            Dictionary of chemical potentials.
        """
        chempots_boundary ={}
        for el,chempot in fixed_chempot.items():
            el_fixed, mu_fixed = el, chempot
        pd = PhaseDiagram(self.get_pd_entries())
        for e in pd.all_entries:
            if e.composition == comp1:
                entry1=e
            if e.composition == comp2:
                entry2=e    
        e1 = pd.get_form_energy(entry1)
        e2 = pd.get_form_energy(entry2)
        
        coeff1 = []
        coeff2 = []
        # order of variables (mu) will follow the order of self.chempots_reference which is alphabetically ordered
        for el in self.chempots_reference:
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
        for el in self.chempots_reference:
            if el != el_fixed:
                chempots_boundary[el] = x[counter]
                counter += 1
                
        return chempots_boundary
                  
                    
        
