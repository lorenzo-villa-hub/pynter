#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:40:26 2020

@author: villa
"""
import pynter
import re
import os
import os.path as op
import json
import pkgutil
import pandas as pd
from monty.json import jsanitize,MontyEncoder, MontyDecoder


def change_file (input_file , output_file=None, back_up_file = True,
                 str_to_modify = {}, lines_to_add=[], check_added_lines = True, lines_to_remove = []):
    
    ''' function to generate a copy of a file changing a string in specific lines
         In case input and output files are the same a copy of the old file is 
         created with name 'old_{input_file}'
    
     input_file: input file to change
     output_file: output file to generate
                  - Default is 'None' - in this case the same is set equal as the input file
                  - In case output_file and input_file have the same name a copy of 
                    the original input file is saved
     str_to_modify: dictionary of strings to be modified - format: {'old string':'new string'}
     lines_to_add: list of lines to be added at the end of the file: format: ['string 1','string 2']
     check_added_lines: - True - if line is already present is not added
                        - False - line is added indipendently of the others
    ''' 
    origin_file = open(input_file,'r')
    lines_list = origin_file.readlines()
    origin_file.close()
    
    if output_file == None or input_file == output_file:
        output_file = input_file
        if back_up_file == True:
            os.rename(input_file, 'old_' + input_file) 
        else:
            os.remove(input_file)
    
    new_file = open(output_file,'w')
           
    for line in lines_list:        
        rem_line = False
        mod_line = False 
        for l in lines_to_remove:
            if l + '\n' == line:
                rem_line = True
        for string in str_to_modify: 
            target_line = re.findall(string, line)                      
            if target_line != [] and rem_line == False:
                new_file.write(line.replace(string,str_to_modify[string]))
                mod_line = True                  
        if mod_line == False and rem_line == False:
             new_file.write(line)   
         
    for l in lines_to_add:
          if check_added_lines:              
              if l + '\n' not in lines_list:            
                  new_file.write(l + '\n')
          else:
              new_file.write(l + '\n')
    
    new_file.close()
    
    return
    

def display_df(df,max_rows=None,max_columns=None,max_colwidth=None):
    """
    Print DataFrame with temporary options

    Parameters
    ----------
    df : (DataFrame)
        DataFrame object.
    max_rows : (int), optional
        Max number of rows to display. The default is None.
    max_columns : (int), optional
        Max number of columns to display. The default is None.
    max_colwidth : (int), optional
        Max column width to display. The default is None.
    """
    with pd.option_context('display.max_rows', max_rows, 'display.max_columns', max_columns,
                           'display.max_colwidth', max_colwidth):  # more options can be specified also
        print(df)
    return


def set_display_df(reset=False,max_rows=None,max_columns=None,max_colwidth=None):
    """
    Set option for DataFrame display.

    Parameters
    ----------
    reset : (bool), optional
        Reset options to default. The default is False.
    max_rows : (int), optional
        Max number of rows to display. The default is None.
    max_columns : (int), optional
        Max number of columns to display. The default is None.
    max_colwidth : (int), optional
        Max column width to display. The default is None.
    """
    if reset:
        pd.reset_option('display.max_rows')
        pd.reset_option('display.max_columns')
        pd.reset_option('display.max_colwidth')
    else:
        pd.set_option('display.max_columns', max_rows)  
        pd.set_option('display.max_rows', max_columns)  
        pd.set_option('display.max_colwidth', max_colwidth)
    return


def explore_packages(path):
    """
    Return package names contained in a path
    """
    names = []
    for sm in pkgutil.walk_packages([path]):
        names.append(sm[1])
        
    return names
    
def explore_pynter_packages():
    path = pynter.__path__[0]
    return explore_packages(path)            


def get_content_from_url(file_url):
    """
    Get content from URL.

    Parameters
    ----------
    file_url: URL to retrieve content from.

    Returns
    -------
    String with content.
    """
    import urllib.request
    with urllib.request.urlopen(file_url) as url:
       content = url.read().decode('utf-8')
    return content
    

def get_object_feature(obj,feature):
    """
    Get value of attribute or method of a generic Object.
    If feature is a single method only the string with the method's name is required.
    If the target feature is stored in a dictionary (or dict of dictionaries), a list of this format needs to be provided:
        ["method_name",key1,key2,...] - This will identify the value of Object.method[key1][key2][...] .

    Parameters
    ----------
    obj : (object)
        Generic object.
    feature : (str or list)
        Method or attribute of class for which the value is needed.
    """
    if isinstance(feature,list):
        method = feature[0]
        try:
            attr = getattr(obj,method) ()
        except:
            attr = getattr(obj,method)                
        for k in feature[1:]:
            try:
                v = attr[k]
            except:
                v = None
            if isinstance(v,dict):
                if feature.index(k) + 1 == len(feature):
                    return v
                else:
                    attr = v
            else:
                return v
            
    else:
        met = feature
        try:
            attr = getattr(obj,met) ()
        except:
            attr = getattr(obj,met)
        return attr


def get_object_from_json(object_class,path_or_string):
    """
    Build class object from json file or string. The class must posses the 'from_dict' method.

    Parameters
    ----------
    object_class : (class)
        Class of the object to decoder.
    path_or_string : (str)
        If an existing path to a file is given the object is constructed reading the json file.
        Otherwise it will be read as a string.

    Returns
    -------
    Decoded object.
    """
    if op.isfile(path_or_string):
        with open(path_or_string) as file:
            d = json.load(file)
    else:
        d = json.loads(path_or_string)

    return object_class.from_dict(d)


def save_object_as_json(object,path,sanitize=False,cls=MontyEncoder):
    """
    Save class object as json string or file. The class must posses the 'as_dict' method.

    Parameters
    ----------
    object: object of a class
    path : (str)
        Path to the destination file.  If None a string is exported.

    Returns
    -------
    d : (str)
        If path is not set a string is returned.
    """
    d = object.as_dict()
    if sanitize:
        d = jsanitize(d)
    if path:
        with open(path,'w') as file:
            json.dump(d,file,cls=cls)
        return
    else:
        return json.dumps(d,cls=cls) 


def select_objects(objects,mode='and',exclude=False,functions=None,**kwargs):
    """
    Filter objects based on different criteria. Returns a list of objects.

    Parameters
    ----------
    objects : (list)
        List of objects.
    mode : (str), optional
        Filtering mode, possibilities are: 'and' and 'or'. The default is 'and'. 
    exclude : (bool), optional
        Exclude the entries satisfying the criteria instead of selecting them. The default is False.
    functions : (list), optional
        Functions containing criteria. The functions must take the object as the argument and
        return a bool. 
    **kwargs : (dict)
        Keys are methods/attributes the objects have. Values contain the criteria. 
        To address more than one condition relative to the same attribute,
        use lists or tuples (e.g. charge=[0,1]).

    Returns
    -------
    output_objects : (list)
        List with selected objects.
    """    
    selected_objects = []
    entered_selection = False
    filtered_objects = objects.copy()
    if functions:
        for func in functions:
            if func:
                entered_selection = True
                if mode=='and':
                    if selected_objects:
                        filtered_objects = selected_objects.copy()
                        selected_objects = []
                for obj in filtered_objects:
                    if func(obj) == True:
                        if obj not in selected_objects:
                            selected_objects.append(obj)
    
    for key in kwargs:
        if mode=='and':
            if selected_objects or entered_selection is True:
                filtered_objects = selected_objects.copy()
                selected_objects = []
        for obj in filtered_objects:
            feature = get_object_feature(obj,key)
            if type(kwargs[key]) in [list,tuple]:
                for value in kwargs[key]:
                    if feature == value:
                        if obj not in selected_objects:
                            selected_objects.append(obj)
            elif feature == kwargs[key]:
                if obj not in selected_objects:
                    selected_objects.append(obj)  
    
    output_objects = []
    for obj in objects:
        if exclude:
            if obj not in selected_objects:
                output_objects.append(obj)
        else:
            if obj in selected_objects:
                output_objects.append(obj)    
    
    return output_objects



def sort_objects(objects,features,reverse=False):
    """
    Sort objects based on a list of features (attributes or methods of the objects, or functions)

    Parameters
    ----------
    objects : (list)
        List of objects to sort.
    features : (list)
        List of features (see get_object_feature).
    reverse : (bool)
        Reverse order.

    Returns
    -------
    (list)
        Sorted objects.
    """
    def criteria(obj):
        criteria = []
        for feature in features:
            if callable(feature):
                criteria.append(feature(obj))
            else:
                criteria.append(get_object_feature(obj,feature))
        return criteria
    sort_function = lambda obj : criteria(obj)
    sorted_objects = sorted(objects,key=sort_function,reverse=reverse)
    
    return sorted_objects


def grep(search_string,file):    
    '''
    Function that emulates "grep" function and returns a list of lines
    that contain target string
    '''         
    lines = []    
    with open (file,'r') as origin_file:
    # searching lines in input file   
        for line in origin_file:
       # emulating 'grep' command
            target_line = re.findall(search_string, line)
            if target_line:
                lines.append(line)
        return lines
         

def grep_list(search_string,target_list):
    """
    Search string in a list of strings. Returns lines that contain the searched string.
    """
    lines=[]
    for l in target_list:
        if search_string in l:
            lines.append(l)
            
    return lines
