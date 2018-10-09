import networkx as nx
import Levenshtein
import numpy as np
from scipy.sparse import coo_matrix
from itertools import combinations_with_replacement, permutations, product
from collections import defaultdict
from functools import partial


def compare_columns(df, func, symmetric=True):
    import pandas as pd

    if symmetric:
        cols = combinations_with_replacement(df.columns, 2)
    else:
        cols = product(df.columns, df.columns)

    df2 = pd.DataFrame([[c1, c2, func(df[c1],df[c2])] for c1,c2 in cols])
    name = df.columns.name
    compare_name = func.__name__
    if name:
        df2.columns = name + '_A', name + '_B', compare_name
    
    if symmetric:
        a,b,c = df2.columns
        dummy = df2[[b,a,c]]
        dummy.columns = a,b,c
        df2 = pd.concat([df2, dummy])[df2.columns]
        return df2.drop_duplicates()
    return df2
     
def compare(df, func=None, symmetric=True, pivot=True):
    """return heatmap looking thing
    """
    if func is None:
        func = jaccard_series
        
    df2 = compare_columns(df, func, symmetric=symmetric)
    if pivot:
        index, col, val = df2.columns
        return df2.pivot(index=index,columns=col,values=val)
    return df2

def in_common_series(a, b):
    """ Assumes inputs are numpy/pandas vectors of counts.
    Counts in a from items that are nonzero in b.
    Not symmetric.
    """
    return a[b>0].sum()

def jaccard(a, b):
    A, B = set(a), set(b)
    return len(A&B) / float(len(A|B))

def jaccard_series(a,b):
    """ Assumes inputs are numpy/pandas vectors of counts.
    """
    A, B = a>0, b>0
    return (A&B).sum() / float((A|B).sum())

def overlap_series(a,b,fraction=True):
    """ Assumes inputs are numpy/pandas vectors of counts.
    """
    A, B = a>0, b>0
    if fraction:
        return (A&B).sum() / float(A.sum())
    return (A&B).sum()

def fold_change_series(a,b):
    """Average log fold change for barcodes present in both series.
    """
    A, B = a>0, b>0
    ix = A | B
    # TODO: fix so input can have zero abundance. 
    return -1 * np.mean(np.abs(np.log10(a[ix] / b[ix])).astype(float))

def pearsonr(a,b):
    import scipy.stats
    return scipy.stats.pearsonr(a,b)[0]

### EDIT DISTANCE

def sparse_edit_distance(sequences):
    """Returns sparse binary matrix of sequence pairs with an edit distance of 1.
    All sequences must be the same length. Matrix is upper triangular.
    """
    overlap = defaultdict(list)
    n = len(sequences[0])
    shape = (len(sequences),)*2

    sequences = [s + s for s in sequences]
    for i in range(n):
        for j, s in enumerate(sequences):
            key = s[i:i+n-1]
            overlap[(i,key)].append(j)

    pairs = [] 
    for vals in overlap.values():
        for entry in combinations(vals, 2):
            pairs.append(entry)
    if pairs:
        data = [True] * len(pairs)
        
        return coo_matrix((data, zip(*pairs)), shape)
    else:
        return coo_matrix(([], ([], [])), shape=shape, dtype=np.bool)

def edit_distance(sequences):
    """Calculate pairwise distances, return in lower triangular matrix.
    """
    distances = np.zeros(2 * [len(sequences)])
    for i, a in enumerate(sequences):
        for j, b in enumerate(sequences[:i]):
            distances[i,j] = Levenshtein.distance(a, b)
    return distances

def make_adjacency(distances, counts, threshold=1):
    """
    Input: pairwise distance (doesn't need to be symmetric?)
    Output: networkx DiGraph
    """
    counts = np.array(counts)
    distances = np.array(distances)
    # find local maxima in adjacency graph
    G = (distances <= threshold) & (distances > 0)

    # scale columns by total reads
    # symmetrize
    G = ((1*G + 1*G.T) > 0).astype(float)
    G = G * counts - counts[:,None]
    G[G < 0] = 0
    return G

def edit_str(a,b):
    arr = []
    op_dict = {'replace': 'R', 'delete': 'D', 'insert': 'I'}
    for op, i, j in Levenshtein.editops(a, b):
        arr += ['%s%d:%s->%s' % (op_dict[op], i, a[i], b[j])]
        # need to track actual insertions/deletions from here on...
        if op in ('insert', 'delete'):
            arr += ['***']
            break
    return '\n'.join(arr)
 
# class UMIGraph(nx.DiGraph):

#     def __init__(self, sequences, counts, threshold=1):
#         self.sequences = sequences
#         self.counts = counts
#         self.distances = edit_distance(sequences)
#         self.G = make_adjacency(self.distances, counts, threshold=threshold)
#         super(UMIGraph, self).__init__(self.G)

#         self.components = list(self.find_components())

#     def find_components(self):
#         components = nx.weakly_connected_component_subgraphs(self)
#         return sorted(list(components), key=lambda x: -1 * len(x.nodes()))

#     def _default_label(self, base, node):
#         edit = edit_str(self.sequences[node], self.sequences[base])
#         return '%s\n%d' % (edit, self.counts[node]) 

#     def _default_base_label(self, base, node):
#         return '%s\n%d' % (self.sequences[base], self.counts[base])
    
#     def label_nodes(self, subgraph, label_fcn=None, base_label_fcn=None):
        
#         label_fcn = label_fcn or self._default_label
#         base_label_fcn = base_label_fcn or self._default_base_label

#         base = nx.topological_sort(subgraph.reverse())[0]
#         labels = {}
#         for node in subgraph.nodes():
#             if node == base:
#                 labels[node] = base_label_fcn(base, node)
#             else:
#                 labels[node] = label_fcn(base, node)

#         return labels

#     def __class__(self):
#         """Needed for self.copy() to work.
#         """
#         return nx.DiGraph().__class__()

#     def draw(self, subgraph, ax=None, labels=None, layout=nx.circular_layout):
#         import matplotlib.pyplot as plt



#         ax = ax or plt.gca()
#         labels = labels or self.label_nodes(subgraph)

#         node_size = np.array(self.counts)[subgraph.nodes()]
#         node_size = (node_size / max(node_size))*2000 + 200
#         print node_size

#         nx.draw_networkx(subgraph, labels=self.label_nodes(subgraph), 
#                          color='w', node_size=node_size, pos=layout(subgraph), 
#                          ax=ax)






