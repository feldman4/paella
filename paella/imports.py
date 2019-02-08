import re
from collections import Counter, OrderedDict, defaultdict
import pandas as pd
from Levenshtein import distance
import numpy as np
from itertools import product, combinations
from natsort import natsorted
import os
from glob import glob
import functools
from functools import partial

import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
sns.set(style='white', font_scale=1.5)

import IPython
from IPython.display import display

IPython.get_ipython().magic('matplotlib inline')
IPython.get_ipython().magic('load_ext autoreload')
IPython.get_ipython().magic('autoreload 2')

from paella import pairwise
import paella.flow
from paella.analysis import *
from paella.utils import and_join, or_join
from tqdm import tqdm_notebook as tqdn

FSC_A  = 'FSC_A'
SSC_A  = 'SSC_A'
FSC_H  = 'FSC_H'
FSC_W  = 'FSC_Width'
SSC_H  = 'SSC_H'
FITC   = 'FITC_A'
ECD    = 'ECD_A'
PE     = 'PE_A'
APC    = 'APC_A'
TIME   = 'Time'
PB450  = 'PB450_H'
KO525  = 'KO525_H'