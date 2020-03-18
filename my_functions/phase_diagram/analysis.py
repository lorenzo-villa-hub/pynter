
import numpy as np
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram, GrandPotentialPhaseDiagram


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
        self._chempots_reference = PDHandler(self._computed_phases).get_chempots_reference()


    @property
    def computed_phases(self):
        return self._computed_phases


    @property
    def chempots_reference(self):
        return self._chempots_reference
    
       
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
                    
            
    def get_chempots_boundary(self,comp1,comp2,fixed_chempot_delta):
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
        fixed_chempot_delta : (Dict)
            Dictionary with fixed Element as key and respective chemical potential as value ({Element:chempot}). The chemical potential
            used here is the one relative to the reference (delta_mu)

        Returns
        -------
        chempots_boundary : (Dict)
            Dictionary of chemical potentials.
        """
        chempots_boundary ={}
        for el,chempot in fixed_chempot_delta.items():
            el_fixed, mu_fixed = el, chempot
        pdhandler = PDHandler(self.computed_phases)
        e1 = pdhandler.get_formation_energy_from_stable_comp(comp1)
        e2 = pdhandler.get_formation_energy_from_stable_comp(comp2)
        
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
                  
    
    def get_composition_boundaries(self,comp,fixed_chempot):
        """
        Get compositions of phases in boundary of stability with a target composition given a fixed chemical potential 
        on one component. Currently only works for 3-component PD (to check). 
        Used Pymatgen GrandPotentialPhaseDiagram class. Careful that the fixed chemical potential has to be the global value,
        not the referenced (mu = mu_ref + delta_mu)

        Parameters
        ----------
        comp : (Pymatgen Composition object)
            Target composition for which you want to get the bounday phases.
        fixed_chempot : (Dict)
            Dictionary with fixed Element as key and respective chemical potential as value ({Element:chempot}).
            The chemical potential here is the global value. 

        Returns
        -------
        comp1,comp2 : (Pymatgen Composition object)
            Compositions of the boundary phases given a fixed chemical potential for one element.
        """
        
        entries = PDHandler(self.computed_phases).pd_entries()
        gpd = GrandPotentialPhaseDiagram(entries, fixed_chempot)
        stable_entries = gpd.stable_entries
        comp_in_stable_entries = False
        for e in stable_entries:
            if e.original_comp == comp:
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
        return comp1,comp2
            
        
        
class PDHandler:
    
    def __init__(self,computed_phases):
        self._computed_phases = {Composition(phase):computed_phases[phase] for phase in computed_phases}

        
    @property
    def computed_phases(self):
        return self._computed_phases

        
    def get_chempots_reference(self):
        """
        Gets elemental reference compounds and respective e.p.a with Pymatgen el_ref attribute in PhaseDiagram class
        the difference with my version is that this considers having more than 1 calculation and takes the minimum e.p.a
        Returns
        -------
        chempot_ref: (Dict)
            Dictionary of elemental chemical potentials 
        """
        
        chempots_ref = {}
        pd = PhaseDiagram(self.pd_entries())
        for el in pd.el_refs:
            chempots_ref[el] = pd.el_refs[el].energy_per_atom
        chempots_ref = {k: v for k, v in sorted(chempots_ref.items(), key=lambda item: item[0])}
        return chempots_ref
       
        
    def get_entries_from_comp(self,comp):
       """
       Get a list of entries corrisponding to the target composition.
       If single_entry is set to True the output a single entry instead of a list
       
       Parameters
       ----------
       comp : (Pymatgen Composition object)
       Returns
       -------
       List of Pymatgen PDEntry objects
       """
       target_entries=[]
       pd = self.phase_diagram()
       for e in pd.all_entries:
           if e.composition.reduced_composition == comp:
               target_entries.append(e)
       if target_entries is not []:
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
        pd = self.phase_diagram()
        form_energies = {}
        for e in self.get_entries_from_comp(comp):
            form_energies[e] = pd.get_form_energy(e)
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
        pd = self.phase_diagram()
        entry = self.get_stable_entry_from_comp(comp)
        return pd.get_form_energy(entry)


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
        pd = self.phase_diagram()
        for e in pd.stable_entries:
            if e.composition == comp:
                target_entry = e
                break
        if target_entry is not None:
            return target_entry
        else:
            raise ValueError('No stable entry has been found for target composition:%s' %comp.reduced_formula)
        
             
    def pd_entries(self):
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
    
    
    def phase_diagram(self):
        """
        Gets Pymatgen PhaseDiagram object 
        Returns
        -------
        Pymatgen PhaseDiagram object
        """
        return PhaseDiagram(self.pd_entries())
    
    