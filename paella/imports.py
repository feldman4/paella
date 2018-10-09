import regex as re
from collections import Counter, OrderedDict, defaultdict
import pandas as pd
from Levenshtein import distance
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from itertools import product, combinations
from natsort import natsorted
import os
from glob import glob
import functools
from functools import partial


sns.set(style='white', font_scale=1.5)

import IPython
from IPython.display import display

IPython.get_ipython().magic('matplotlib inline')
IPython.get_ipython().magic('load_ext autoreload')
IPython.get_ipython().magic('autoreload 2')

from paella import pairwise

from paella.analysis import *

from tqdm import tqdm

FSC_A  = 'FSC_A'
SSC_A  = 'SSC_A'
FSC_H  = 'FSC_H'
FSC_W  = 'FSC_Width'
SSC_H  = 'SSC_H'
FITC   = 'FITC_A'
ECD    = 'ECD_A'
TIME   = 'Time'

### UTILS


def bin_join(xs, symbol):
    symbol = ' ' + symbol + ' ' 
    return symbol.join('(%s)' % x for x in xs)
        
or_join  = functools.partial(bin_join, symbol='|')
and_join = functools.partial(bin_join, symbol='&')
