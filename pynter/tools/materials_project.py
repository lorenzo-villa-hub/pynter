#!/usr/bin/env python

from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.inputs import Poscar
import argparse as ap
from pynter.__init__ import load_config

API_KEY = load_config()['API_KEY']

with MPRester(API_KEY) as mpr:

    class MPDatabase:
        
        def __init__(self,mp_id=None):
            """
            Class to retrieve data from Materials Project database

            Parameters
            ----------
            mp_id : (str), optional
                Materials-ID. The default is None.
            """
            
            self.mp_id = mp_id if mp_id else None
            self.api_key = API_KEY
        
        @property
        def mp_rester(self):
            return mpr
            
            
        def args(self):
            """
            Get arguments for use from terminal

            Returns
            -------
            args :
                ArgumentParser object.
            """
            
            parser = ap.ArgumentParser()
            
            parser.add_argument('-id','--mp-id',help='Materials project ID',required=True,type=str,metavar='',dest='mp_id')
            parser.add_argument('-s','--structure',help='print Structure',action='store_true',required=False,default=False,dest='print_structure')
            parser.add_argument('-P','--poscar',help='create POSCAR',action='store_true',required=False,default=False,dest='write_poscar')
            parser.add_argument('-c','--conv-cell',help='Get conventional unit cell, Default is primitive',action='store_true',required=False,default=False,dest='get_conventional')
            parser.add_argument('-f','--final-structure',help='Get final structure, Default is False',action='store_true',required=False,default=False,dest='get_final')     
            parser.add_argument('-n','--filename',help='Filename for POSCAR, default is "POSCAR"',required=False,default='POSCAR',type=str,metavar='',dest='filename')
            
            args =  parser.parse_args()
            return args
            
            
        def get_data(self,data_type='vasp'):    
            """
            Flexible method to get any data using the Materials Project REST
            interface. Generally used by other methods for more specific queries.
    
            Format of REST return is *always* a list of dict (regardless of the
            number of pieces of data returned. The general format is as follows:
    
            [{"material_id": material_id, "property_name" : value}, ...]
    
            This is generally a call to
            https://www.materialsproject.org/rest/v2/materials/vasp/<prop>.
            See https://github.com/materialsproject/mapidoc for details.
    
            Parameters
            ----------
                data_type (str): Type of data to return. Currently can either be
                    "vasp" or "exp".
            Returns:
                List of dictionaries with available MP data
            """
            data = mpr.get_data(self.mp_id,data_type=data_type) 
            
            return data


        def get_structure(self,final=False,conventional_unit_cell=False):
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
            return mpr.get_structure_by_material_id(self.mp_id,final=final,conventional_unit_cell=conventional_unit_cell).get_sorted_structure()  

        
if __name__ == '__main__':
    
    mp = MPDatabase()
    args = mp.args()
    mp.mp_id = args.mp_id    
    structure = mp.get_structure(final=args.get_final,conventional_unit_cell=args.get_conventional)
    
    if args.print_structure:
        print(structure)
    if args.write_poscar:
        Poscar(structure).write_file(args.filename)
    
    
    
    
    
    
    
