#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GIFT_lib

GIFT_lib is a python package to run the GIFT algorithm.
https://github.com/songpeng/GIFT
TODO: list of the functions and classes in this module.
"""

# import should from general to specific way.
import argparse
import numpy
import pubchemfp # TODO: get compounds' pubchem fingerprints
import rdkitfp # TODO: get compounds' rdkit based fingerprints
import chemfp.bitops as chembit # fast calculate the jaccard disctance

__author__ = "Songpeng Zu"
__copyright__ = "Copyright 2015, Bioinformatics Lab, Department of Automation, Tsinghua University"
__version__ = "1.1"
__maintainer__ = "Songpeng Zu"
__email__ = "zusongpeng@gmail.com"
__status__ = "Development"

class gift_em:
    """Class to represent the elements needed in the algorithm.
    """
    def __init__(self,snum,dnum):
        pass
