#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 12:06:45 2023

@author: villa
"""

def get_dict_from_line_string(dict_string,convert_to_float=True):
    items = dict_string.split()
    if len(items) % 2 != 0:
        raise ValueError("Invalid dictionary string")

    dictionary = {}
    for i in range(0, len(items), 2):
        key = items[i]
        value = items[i + 1]
        dictionary[key] = float(value) if convert_to_float else value

    return dictionary

def round_floats(obj):
    if isinstance(obj, float):
        return round(obj, 2)
    elif isinstance(obj, dict):
        return {key: round_floats(value) for key, value in obj.items()}
    elif isinstance(obj, list) or isinstance(obj, tuple):
        return [round_floats(element) for element in obj]
    else:
        return obj
