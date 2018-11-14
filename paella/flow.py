import pandas as pd
import numpy as np
import os
import re 

import FlowCytometryTools
import seaborn as sns
import matplotlib.pyplot as plt
from paella.imports import and_join, or_join

def load_fcs(filename):
    fcm = FlowCytometryTools.FCMeasurement(ID='', datafile=filename)
    df = fcm.data
    df = df.rename(columns=lambda x: x.replace('-', '_'))
    return df.assign(file=filename)


def add_sample_info(df):
    """Get sample info from Cytoflex naming scheme.
    """
    pat = '(\d\d)-(?i)Tube-(.*).fcs'
    plate, well = re.findall(pat, df['file'].iloc[0])[0]
    row, col = well[0], int(well[1:])
    return df.assign(plate=plate, well=well, row=row, col=col)


def transform_columns(df, columns, transform):
    df = df.copy()
    df[columns] = transform(df[columns])
    return df


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


def add_gates(df, gates):
    df = df.copy() # copy if we are changing original
    for gate_name in gates:
        df[gate_name] = df.eval(gates[gate_name])
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
    