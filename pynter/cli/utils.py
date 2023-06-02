#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 12:06:45 2023

@author: villa
"""

def get_dict_from_line_string(dict_string,convert_to_float=True):
    """
    Convert a string in the format "key1 value1 key2 value2 ..." to a dictionary.

    Parameters
    ----------
    dict_string : (str)
        String in the format "key1 value1 key2 value2 ...".
    convert_to_float : (bool), optional
        Convert values to float. The default is True.

    Returns
    -------
    dictionary : (dict)
    """
    items = dict_string.split()
    if len(items) % 2 != 0:
        raise ValueError("Invalid dictionary string")

    dictionary = {}
    for i in range(0, len(items), 2):
        key = items[i]
        value = items[i + 1]
        dictionary[key] = float(value) if convert_to_float else value

    return dictionary


def round_floats(obj,decimals=2):
    """
    Round floats in an object (float, dict, tuple or list)

    Parameters
    ----------
    obj : 
        Object.
    decimals : (int)
        Number of decimals to round to.

    Returns
    -------
    Object with rounded floats.

    """
    if isinstance(obj, float):
        return round(obj, decimals)
    elif isinstance(obj, dict):
        return {key: round_floats(value) for key, value in obj.items()}
    elif isinstance(obj, list) or isinstance(obj, tuple):
        return [round_floats(element) for element in obj]
    else:
        return obj
