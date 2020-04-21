#!/nfshome/villa/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 18:04:58 2019

@author: villa
"""

import sys

def grep(search_string,file):
    
    '''Function that emulates "grep" function and returns a list of lines
       that contain target string'''     
    
    import re
    
    lines = []    
    with open (file,'r') as origin_file:
    # searching lines in input file   
        for line in origin_file:
       # emulating 'grep' command
            target_line = re.findall(search_string, line)
            if target_line:
                lines.append(line)
        return lines
         
            
if __name__ == '__main__':
    
    search_string = sys.argv[1]
    file = sys.argv[2]
    
    lines = grep(search_string,file)
    for l in lines:
        print(l)
        
