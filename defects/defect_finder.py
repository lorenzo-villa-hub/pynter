#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 19:08:20 2020

@author: lorenzo
"""
import numpy as np
from pymatgen.core.periodic_table import Element

def defect_finder(structure_defect,structure_bulk): # change with coordinate comparison > more accurate

    df = structure_defect
    bk = structure_bulk
    
    if len(df) > len(bk):
        defect_type = 'Interstitial'
        for s in df:
            if s not in bk:
                defect_site = s
                defect_specie = s.specie
        
    elif len(df) < len(bk):
        defect_type = 'Vacancy'
        for s in bk:
            if s not in df:
                defect_site = s
                defect_specie = s.specie

    elif len(df) == len(bk):
        comp_df = structure_defect.composition
        comp_bk = structure_bulk.composition
        if comp_df != comp_bk:
            defect_type = 'Substitution'
            for el in comp_df:
                if el not in comp_bk:
                    sub_element = el
            for s in df:
                if s.specie == sub_element:
                    defect_site = s
                    defect_specie = s.specie

    return defect_type, defect_site                    
                          

    # if len(df) > len(bk):
    #     defect_type = 'Interstitial'
    #     for sd in df:
    #         site_found = False
    #         df_coords = np.around(sd.frac_coords,decimals=3)
    #         for sb in bk:
    #             bk_coords = np.around(sb.frac_coords,decimals=3)
    #             if df_coords != bk_coords and sb.specie == sb.specie:
    #                 site_found = True
    #                 break
    #         if site_found:
    #             defect_site = sd
    #             return defect_site,defect_type    
        
    # elif len(df) < len(bk):
    #     defect_type = 'Vacancy'
    #     for sb in bk:
    #         site_found = False
    #         bk_coords = np.around(np.asarray(sb.frac_coords,dtype = np.float64),decimals=4)
    #         for sd in df:
    #             df_coords = np.around(np.asarray(sd.frac_coords,dtype = np.float64),decimals=4)
    #             if np.array_equiv(bk_coords,df_coords) and sb.specie == sd.specie:
    #                 site_found = True
    #                 break
    #         if site_found:
    #             defect_site = sb
    #             return defect_site,defect_type 
            
            
    # elif len(df) == len(bk):
    #     comp_df = structure_defect.composition
    #     comp_bk = structure_bulk.composition
    #     if comp_df != comp_bk:
    #         defect_type = 'Substitution'
    #         for el in comp_df:
    #             if el not in comp_bk:
    #                 sub_element = el
    #         for s in df:
    #             if s.specie == sub_element:
    #                 defect_site = s
    #                 defect_specie = s.specie
                    
                           
    # return defect_type, defect_site