import pandas as pd
import numpy as np
import os
import re 

import FlowCytometryTools
import seaborn as sns
import matplotlib.pyplot as plt
from paella.utils import and_join, or_join
import collections

from bs4 import BeautifulSoup

def load_fcs(filename):
    fcm = FlowCytometryTools.FCMeasurement(ID='', datafile=filename)
    df = fcm.data
    df = df.rename(columns=lambda x: x.replace('-', '_'))
    return df.assign(file=filename)


def add_sample_info(df):
    """Get sample info from Cytoflex naming scheme. Slow in spite of caching.
    """
    @Memoized
    def get_info(f):
        pat = '(\d\d)-Tube-(.*).fcs'
        plate, well = re.findall(pat, f)[0]
        row, col = well[0], int(well[1:])
        return plate, well, row, col

    plate, well, row, col = zip(*[get_info(f) for f in df['file']])
    return df.assign(plate=plate, well=well, row=row, col=col)


def parse_sample_info(f):
    pat = '(\d\d)-Tube-(.*).fcs'
    plate, well = re.findall(pat, f)[0]
    row, col = well[0], int(well[1:])
    return dict(plate=plate, well=well, row=row, col=col)



def transform_columns(df, columns, transform):
    # why won't the nice way work??
    df = df.copy()
    for c in columns:
        df[c] = transform(df[c])
    return df
    # return df.assign(**{c: lambda x: transform(x[c]) for c in columns})
    

def plot_flow(df, x, y, ax=None, color=None):
    """throw away color (from seaborn)
    """
    colors = bilinear_interpolate(df[x], df[y])
    if ax is None:
        ax = plt.gca()
    ax.scatter(df[x], df[y], c=colors, s=5, lw=0, cmap='jet')
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    return ax


def plot_flow_sns(**kwargs):
    kwargs['df'] = kwargs.pop('data')
    kwargs.pop('color')
    return plot_flow(**kwargs)


def multi_and_eval(df, gate):
    n = 5
    lines = gate.split('&')
    subgates = [and_join(lines[n*i:n*(i+1)]) for i in range(int(len(lines)/n) + 1)]
    arr = [df.eval(s) for s in subgates if s]
    if len(arr) == 1:
        return arr[0]
    else:
        for x in arr[1:]:
            arr[0] = arr[0] & x
        return arr[0]


def add_gates(df, gates):
    df = df.copy() # copy if we are changing original
    for gate_name in gates:
        df[gate_name] = multi_and_eval(df, gates[gate_name])
    return df


def plot_with_dims(df, x, y, limits=None, ax=None):
    """Plot and apply axis formatting.
    """
    ax = plot_flow(df, x, y, ax=ax)
    if limits:
        for dimension in limits:
            if x == dimension:
                ax.set_xlim(limits[dimension])
            if y == dimension:
                ax.set_ylim(limits[dimension])
    return ax


def bilinear_interpolate(x, y, bins=None):
    """
    FROM FCM PACKAGE

    Returns interpolated density values on points (x, y).
    
    Ref: http://en.wikipedia.org/wiki/Bilinear_interpolation.
    """
    if bins is None:
        bins = int(np.sqrt(len(x)))

    z, unused_xedge, unused_yedge = np.histogram2d(y, x, bins=[bins, bins],
                                        range=[(np.min(y), np.max(y)),
                                               (np.min(x), np.max(x))]
                                        )
    xfrac, xint = np.modf((x - np.min(x)) /
                             (np.max(x) - np.min(x)) * (bins - 1))
    yfrac, yint = np.modf((y - np.min(y)) /
                             (np.max(y) - np.min(y)) * (bins - 1))

    xint = xint.astype('i')
    yint = yint.astype('i')

    z1 = np.zeros(np.array(z.shape) + 1)
    z1[:-1, :-1] = z

    # values at corners of square for interpolation
    q11 = z1[yint, xint]
    q12 = z1[yint, xint + 1]
    q21 = z1[yint + 1, xint]
    q22 = z1[yint + 1, xint + 1]

    return q11 * (1 - xfrac) * (1 - yfrac) + q21 * (1 - xfrac) * (yfrac) + \
        q12 * (xfrac) * (1 - yfrac) + q22 * (xfrac) * (yfrac)
    

class Memoized(object):
    """
    from PythonDecoratorLibrary
    Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value
    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__
    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        return functools.partial(self.__call__, obj) 


## gate extraction

def mx_plus_b(x0, x1, y0, y1):
    m = (y1 - y0) / (x1 - x0)
    b = y0 - m * x0
    return m, b

def make_gate(x0, x1, m, b, kdims):
    comp = '>' if x1 > x0 else '<'
    
    gate = \
    ('y | m*x + b'
     .replace('m', str(m))
     .replace('b', str(b))
     .replace('x', kdims[0])
     .replace('y', kdims[1])
     .replace('|', comp)
    )
    return '(' + gate + ')'

def gates_from_poly(poly, kdims):
    xs = poly['x'][0]
    ys = poly['y'][0]
    return gates_from_xy(xs, ys, kdims)

def gates_from_xy(xs, ys, kdims):
    n = len(xs)
    gates = []
    #print("xs = {}".format(xs))
    for i in range(n):
        x0, x1 = xs[i], xs[(i + 1) % n]
        y0, y1 = ys[i], ys[(i + 1) % n]        
        m, b = mx_plus_b(x0, x1, y0, y1)
        gate = make_gate(x0, x1, m, b, kdims)
        gates += [gate]
    final_gate = ' & '.join(gates)
    return final_gate

def extract_polygon_gates(xml_file):
    with open(xml_file, 'r') as fh:
        xml = fh.read()
    soup = BeautifulSoup(xml, 'xml')
    gates = soup.find_all('PolygonGate')
    
    
    return gates
    
def tag_to_coords(gate_tag):
    arr = []
    for vertex in gate_tag.find_all('vertex'):
        x, y = vertex.find_all('coordinate')
        x = float(x['data-type:value'])
        y = float(y['data-type:value'])
        arr += [(x, y)]

    x, y = gate_tag.find_all('fcs-dimension')
    x = x['data-type:name']
    y = y['data-type:name']
    kdims = [x, y]
    return arr, kdims

def coords_to_gate(xy, kdims):
    xs, ys = zip(*xy)
    kdims = map(lambda x: x.replace('-', '_'), kdims)
    
    return gates_from_xy(xs, ys, kdims)

def log_transform(x, clip=False):
    if clip:
        return np.log10(1. + np.clip(x, 0, None))
    else:
        return np.log10(1 + np.abs(x))

def compress_below(x, cutoff, factor):
    """Linear stand-in for proper transform
    """
    mask = x<cutoff
    x = x.copy()
    x[mask] = (x[mask] - cutoff) / factor + cutoff
    return x

def xml_to_pandas_gates(xml_file, transform=log_transform):
    gate_tags = extract_polygon_gates(xml_file)
    coords = map(tag_to_coords, gate_tags)
    # print("=" * 80)
    # print("COO")
    # print(coords)

    gates = []
    for xy, kdims in coords:
        # guess direction
        x0, y0 = xy[0] 
        x1, y1 = xy[1]
        x2, y2 = xy[2]
        #print(x0, x1, x2)

        dx0, dy0 = x1 - x0, y1 - y0
        dx1, dy1 = x2 - x1, y2 - y1
        cross_product_z = dx0 * dy1 - dx1 * dy0
        #print("=" * 80)
        #print(xy)
        #print("=" * 80)

        if cross_product_z < 0:
            xy = xy[::-1]

        if transform:
            xy = transform(np.array(xy))
        gates += [coords_to_gate(xy, kdims)]
    return gates

