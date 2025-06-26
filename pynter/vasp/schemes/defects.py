#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 13:57:21 2025

@author: villa
"""
import os.path as op

from pynter.defects.defects import create_vacancies, create_substitutions
from pynter.defects.structure import defect_finder
from pynter.vasp.schemes.core import InputSets
from pynter.vasp.schemes.relaxation import RelaxationSchemes
from pynter.vasp.schemes.variations import VariationsSchemes


class DefectSchemes(InputSets):
            
    
    def defects_pbe_relaxation(self,defects_with_charges,locpot=False,rel_scheme='default'):
        """
        Generate jobs for default defect calculation scheme with PBE.

        Parameters
        ----------
        defects_with_charges : (list)
            List of tuples with Defect object and relative charge state list, for example [(Defect , [q1,q2,q3,...])].
        locpot : (bool), optional
            Add 'LVTOT=True' in INCAR. The default is False.
        rel_scheme : (str), optional
            Relaxation scheme to use with PBE. 'default for pbe_relaxation and 
            'gamma' pbe_relaxation_gamma. The default is 'default'.

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        jobs = []

        for df, charge_states in defects_with_charges:
            defect_name = df.name.replace('(','_').replace(')', '') #replace label parenthesis
            path = op.join(self.path,defect_name)
            schemes_q = VariationsSchemes(path=path,structure=df.defect_structure,incar_settings=self.incar_settings,
                                job_settings=self.job_settings,name=self.name+'_'+defect_name,add_parent_folder=False)
            charge_jobs = schemes_q.charge_states(charge_states,locpot)
            for jq in charge_jobs:
                schemes_rel = RelaxationSchemes(path=jq.path,vaspinput=jq.inputs,job_settings=jq.job_settings,name=jq.name)
                if rel_scheme == 'default':
                    rel_jobs = schemes_rel.pbe_relaxation()
                elif rel_scheme == 'short':
                    rel_jobs = schemes_rel.pbe_ionic_rel(add_to_job_name=None,add_to_job_path=None)
                elif rel_scheme == 'gamma':
                    rel_jobs = schemes_rel.pbe_relaxation_gamma()
                for jrel in rel_jobs:
                        jobs.append(jrel)
        return jobs


    def defects_pbe_relaxation_from_structures(self,structures_with_charges,tol=1e-03,locpot=False,rel_scheme='default'):
        """
        Generate jobs for default defect calculation scheme with PBE starting from defect structures.
        defect_finder is used to create defect objects, check that the algorithm wirks as intended.

        Parameters
        ----------
        structures_with_charges : (list)
            List of tuples with Structure object and relative charge state list, for example [(Structure , [q1,q2,q3,...])].
        tol: (float)
            Tolerance for defect_finder.
        locpot : (bool)
            Add 'LVTOT=True' in INCAR. The default is False.
        rel_scheme : (str)
            Relaxation scheme to use with PBE. 'default' for pbe_relaxation and 
            'gamma' pbe_relaxation_gamma. The default is 'default'.

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        defects_with_charges = []
        for structure_defect, charge_states in structures_with_charges:
            defect = defect_finder(structure_defect, self.structure, tol=tol)
            defects_with_charges.append( (defect,charge_states) )
        return self.defects_pbe_relaxation(defects_with_charges,locpot=locpot,rel_scheme=rel_scheme)


    def interstitials_pbe_relaxation(self):
        """
        Since the interstitials structures neeed to be checked and selected, calculations should be set up manually 
        with defects_pbe_relaxation. If the function needs to be implemented, set up a criterion for the choice of 
        interstitials sites.
        """
        pass

            
    def substitutions_pbe_relaxation(self,elements_to_replace_with_charges,supercell_size=None,locpot=False,rel_scheme='default'):
        """
        Generate jobs for default substitutions calculation scheme with PBE.

        Parameters
        ----------
        elements_to_replace_with_charges : (dict)
            Dictionary which indicates the substitution types and the relative charge states list.
            Needs to be in this format:
                {'<new_el>-on-<old_el>':[q0,q1,q2]}, where <new_el> and <old_el> are the symbols of the element to be added and
                the element to be substituted with, respectively. [q0,q1,q2] is the list of corrisponding charge states for this substitution.
            The format of the substitution string is the same that is generated from the function "create_substitutiion_structures".
        supercell_size : (int or numpy array), optional
            Input for the make_supercell function of the Structure class. The default is 1 (the structure is not modified).
        locpot : (bool), optional
            Add 'LVTOT=True' in INCAR. The default is False.
        rel_scheme : (str), optional
            Relaxation scheme to use with PBE. 'default for pbe_relaxation and 
            'gamma' pbe_relaxation_gamma. The default is 'default'.

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        defects_with_charges = []
        for key,charge_states in elements_to_replace_with_charges.items():
            elements_to_replace = {key.split('-on-')[1] : key.split('-on-')[0]}
            substitution = create_substitutions(self.structure,elements_to_replace=elements_to_replace,
                                                supercell_size=supercell_size)[0]
            defects_with_charges.append((substitution,charge_states))
        jobs = self.defects_pbe_relaxation(defects_with_charges,locpot=locpot,rel_scheme=rel_scheme)
        return jobs        
        

    def vacancies_pbe_relaxation(self,elements_with_charges,supercell_size=None,locpot=False,rel_scheme='default'):
        """
        Generate jobs for default vacancies calculation scheme with PBE.

        Parameters
        ----------
        elements_with_charges : (dict)
            Dictionary with element symbols as keys and relative charge states (int) lists as values ({'el':[-1,0,1]}).
        supercell_size : (int or numpy array), optional
            Input for the make_supercell function of the Structure class.
            If None the input structure is not modified. The default is None.
        locpot : (bool), optional
            Add 'LVTOT=True' in INCAR. The default is False.
        rel_scheme : (str), optional
            Relaxation scheme to use with PBE. 'default for pbe_relaxation and 
            'gamma' pbe_relaxation_gamma. The default is 'default'.

        Returns
        -------
        jobs : (list)
            List of VaspJob objects.
        """
        defects_with_charges = []
        for el,charge_states in elements_with_charges.items():
            vacancy = create_vacancies(structure=self.structure,elements=[el],supercell_size=supercell_size)[0]
            defects_with_charges.append((vacancy,charge_states))
        jobs = self.defects_pbe_relaxation(defects_with_charges,locpot=locpot,rel_scheme=rel_scheme)
        return jobs