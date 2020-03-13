#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 15:43:14 2020

@author: villa
"""

def kpoints_file_for_HSE(k_points, k_path, npoints=10, input_file='IBZKPT' , output_file='KPOINTS'):
    ''' Creates KPOINTS file in the format for HSE bandstructure calculations on VASP
        HSE bandstructure calculations on VASP need a KPOINTS file in the format of a IBZKPT file
        The kpoints of the previous SCF calculation need to be preserved, and the kpoints for 
        bandstructure calculations need to be added with weight zero.
        This function adds to the IBZKPT file of the previous SCF calc the interpolation between
        the given kpoints for the kpath
        Args:
            k_points: Dictionary with as key the label of k point and as value the coordinates as a numpy array
            k_path: List of lists of labels of k points that give the order of the path
                    Every list gives which point interpolation is performed, between different
                    lists there is discontinuity (jump in the path)
            npoints: Default(npoints=10)
                    number of k-points interpolated between each special point of the path
            input_file: IBZKPT format file to add kpoints to
            output_file: Output KPOINTS file for HSE calc
    '''
    import numpy as np

    k_list = []    
    # getting coordinates of the path
    for list in k_path:
        k_list_single = []
        for item in list:
            k_list_single.append(k_points[item])
        k_list.append(k_list_single)
    
    k_inter = []   
    # interpolation in between each point of the path
    for l in k_list:
        for p in range(0,len(l)-1):
            k_inter.append(l[p])
            # calculating step for signle point - all coordinates are treated at the same time
            step = (l[p + 1] - l[p])/npoints
            for i in range(1,npoints):
                # adding the i-th step
                point = l[p] + step*i
                k_inter.append(point)
      #      if l == k_list[-1]:
            if p == len(l)-2:
                    k_inter.append(l[p+1])
            # resetting step 
            step = 0
    # getting number of added points
    total_new_kpoints = len(k_inter)    
    
    # creating new KPOINTSfile
    k_file = open(output_file,'w+')
    
    # getting last line of IBZKPT file
    with open(input_file) as origin_file:
        all_lines = origin_file.readlines()
    last_line = all_lines[-1]
    
    # writing new KPOINTS file
    # copying IBZKPT file exept for 2nd line where the total number of kpoints is updated
    # adding after last line of IBZKPT the new kpoints with weight 0
    n_line = 0
    with open(input_file) as origin_file:
        for line in origin_file:
            n_line += 1
            # if line is the 2nd line update total number of k-points and write
            if n_line == 2:
                total_old_kpoints = int(line)
                total_kpoints = total_old_kpoints + total_new_kpoints
                k_file.write('%i \n' %(total_kpoints))
            else:
                k_file.write(line)
            # if line is the last line write new kpoints with weight 0    
            if line == last_line:
                for p in range(0,len(k_inter)):
                    # if point is a special point add label
                    k_labeled = False
                    for key in k_points:
                            if np.all(k_inter[p] == k_points[key]) and k_labeled == False:
                                k_file.write('%f %f %f 0 %s \n' %(k_inter[p][0] , k_inter[p][1] , k_inter[p][2] , key))
                                k_labeled = True
                    if k_labeled == False:            
                        k_file.write('%f %f %f 0 \n' %(k_inter[p][0] , k_inter[p][1] , k_inter[p][2]))
                
    k_file.close()
    
    return

