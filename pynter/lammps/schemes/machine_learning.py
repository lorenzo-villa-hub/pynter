#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  8 17:48:31 2025

@author: villa
"""

from pymatgen.io.lammps.inputs import LammpsInputFile

from pynter.lammps.schemes.core import InputSets


class PACEActiveLearningSchemes(InputSets):
    
    
    def cook_and_quench(self,
                        potential_filename,
                        asi_filename,
                        elements_string=None,
                        add_to_job_name=None,
                        add_to_job_path=None):
        
        elements_string = elements_string or ' '.join(str(el) for el in self.structure.composition)
        input_string = f"""
        units       metal
        boundary    p p p
        atom_style  atomic
        box tilt large        
        read_data    {self.data_filename}
    
        pair_style  pace/extrapolation
        pair_coeff  * * {potential_filename} {asi_filename} {elements_string}
        
        fix pace_gamma all pair 10 pace/extrapolation gamma 1   
        compute max_pace_gamma all reduce max f_pace_gamma      
                     
        variable dump_skip equal "c_max_pace_gamma < 5"         
        dump pace_dump all custom 20 extrapolative_structures.dump id type x y z f_pace_gamma
        dump_modify pace_dump skip v_dump_skip
        
        variable max_pace_gamma equal c_max_pace_gamma         
        fix extreme_extrapolation all halt 10 v_max_pace_gamma > 25
        
        
        thermo 500    
        timestep 0.001
        
        
        velocity all create 100.0 42
        fix 1 all nvt temp 50 300 0.1 
        run 10000
        unfix 1
        
        
        fix 1 all npt temp 300 300 0.1 aniso 0 10000 1
        run 20000
        unfix 1
        
        fix 1 all npt temp 300 300 0.1 aniso 10000 0 1
        run 20000
        unfix 1
        
        
        fix 1 all nvt temp 300 2000 0.1 
        run 68000
        unfix 1
        
        fix 1 all nvt temp 2000 300 0.1 
        run 68000
        unfix 1
                                
        fix 1 all box/relax aniso 0.0
        minimize  0 1.0e-8 10000 1000000
        unfix 1
        
        write_data out.data
        """
        self.lammps_input = LammpsInputFile.from_str(input_string)
        return self.get_lammpsjob(add_to_job_name=add_to_job_name,add_to_job_path=add_to_job_path)
            
