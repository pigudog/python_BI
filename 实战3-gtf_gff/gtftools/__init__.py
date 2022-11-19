# -*- coding: utf-8 -*-
"""
日期：19/11/2022 13:42
作者： PigUDog
"""
import os.path
import sys
from os import path
d = path.dirname(__file__)
parent_path = os.path.dirname(d)
sys.path.append(parent_path)
from gtftools._version import version as __version__
from gtftools._version import author as __author__
from gtftools._version import email as __email__
from gtftools._error import ParsingError
from gtftools._parser import Gtf, Gff3