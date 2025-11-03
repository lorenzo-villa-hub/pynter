#!/usr/bin/env python

#from pymatgen.ext.matproj import MPRester
from mp_api.client import MPRester


class  MPDatabase:
    
    def __init__(self,mp_id=None,API_KEY=None):
        """
        Class to retrieve data from Materials Project database

        Parameters
        ----------
        mp_id : (str), optional
            Materials-ID. The default is None.
        """
        
        self.mp_id = mp_id if mp_id else None
        self.API_KEY = API_KEY
    
    @property
    def mp_rester(self):
        return MPRester(self.API_KEY)
        
        
    def get_entries(self,
                    chemsys_formula_mpids,
                    compatible_only=True,
                    property_data=None,
                    conventional_unit_cell=False):
        """
        Get a list of ComputedEntries or ComputedStructureEntries corresponding
        to a chemical system, formula, or materials_id or full criteria.

        Parameters
        ----------
            chemsys_formula_id_criteria (str/dict): A chemical system
                (e.g., Li-Fe-O), or formula (e.g., Fe2O3) or materials_id
                (e.g., mp-1234) or full Mongo-style dict criteria.
            compatible_only (bool): Whether to return only "compatible"
                entries. Compatible entries are entries that have been
                processed using the MaterialsProject2020Compatibility class,
                which performs adjustments to allow mixing of GGA and GGA+U
                calculations for more accurate phase diagrams and reaction
                energies.
            property_data (list): Specify additional properties to include in
                entry.data. If None, no data. Should be a subset of
                supported_properties.
            conventional_unit_cell (bool): Whether to get the standard
                conventional unit cell.

        Returns:
            List of ComputedEntry or ComputedStructureEntry objects.
        """
        with MPRester(self.API_KEY) as mpr:
            entries = mpr.get_entries(
                                chemsys_formula_mpids=chemsys_formula_mpids,
                                compatible_only=compatible_only,
                                property_data=property_data,
                                conventional_unit_cell=conventional_unit_cell)
        return entries
    
    
    def get_entries_from_compositions(self,
                                    compositions,
                                    compatible_only=True,
                                    property_data=None,
                                    conventional_unit_cell=False):
        """
        Get a dictionary with compositions (strings) as keys and a list of ComputedEntries 
        or ComputedStructureEntries as values.

        Parameters
        ----------
            compositions (list): List of strings with compositions.
            stable_only (bool): Get phase with lowest E above hull.
            compatible_only (bool): 
                Whether to return only "compatible"
                entries. Compatible entries are entries that have been
                processed using the MaterialsProject2020Compatibility class,
                which performs adjustments to allow mixing of GGA and GGA+U
                calculations for more accurate phase diagrams and reaction
                energies.
            property_data (list):
                Specify additional properties to include in
                entry.data. If None, no data. Should be a subset of
                supported_properties.
            conventional_unit_cell (bool):
                Whether to get the standard conventional unit cell

        Returns:
            List of ComputedEntry or ComputedStructureEntry objects.
        """
        entries_dict = {}
        for comp in compositions:
            entries = self.get_entries(
                                    chemsys_formula_mpids=comp,
                                    compatible_only=compatible_only,
                                    property_data=property_data,
                                    conventional_unit_cell=conventional_unit_cell)
            entries_dict[comp] = entries
        
        return entries_dict
    

    def get_structure(self,final=True,conventional_unit_cell=False):
        """
        Get a Structure corresponding to a material_id.

        Parameters
        ----------
            material_id (str): Materials Project material_id (a string,
                e.g., mp-1234).
            final (bool): Whether to get the final structure, or the initial
                (pre-relaxation) structure. Defaults to True.
            conventional_unit_cell (bool): Whether to get the standard
                conventional unit cell

        Returns:
            Structure object.
        """
        with MPRester(self.API_KEY) as mpr:
            structure = mpr.get_structure_by_material_id(self.mp_id,final=final,conventional_unit_cell=conventional_unit_cell)
        return structure

    
    
    
    
    
    
    
