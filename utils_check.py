#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 12:51:32 2017

@author: paulolomeu
"""

from collections import Mapping, Container
from sys import getsizeof

def deep_getsizeof(o, ids): #Find the memory footprint of a Python object

    #:param o: the object
    #:param ids:
    #:return:

    d = deep_getsizeof
    
    if id(o) in ids:
        return 0
    
    r = getsizeof(o)
    ids.add(id(o))
    
    if isinstance(o, str):
        return r
    
    if isinstance(o, Mapping):
        return r + sum(d(k, ids) + d(v, ids) for k, v in o.iteritems())
    
    if isinstance(o, Container):
        return r + sum(d(x, ids) for x in o)
    
    return r