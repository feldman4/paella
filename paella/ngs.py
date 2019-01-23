import numpy as np
import pandas as pd
from natsort import natsorted

def load_info(f):
    """Load sample info from deconvolution tab.
    """
    get_columns = lambda x: x[[c for c in x.columns if 'Unnamed' not in c]]
    df_info = (pd.read_csv(f)
     .pipe(get_columns) # get just the named columns
     .dropna(subset=['ID']).set_index('ID') # set library ID as index
     .rename(columns=lambda x: x.replace(' ', '_').lower())
     )

    df_info['name'] = (df_info['cell_type'] + '_' + 
                   df_info['sample_info'] + '_' +
                   df_info.index)

    return df_info

# bad_barcodes = pd.read_csv('bad_barcodes.csv', header=None)[0].pipe(list)

def load_hist(f, bad_seqs=None, threshold=3):
    if bad_seqs is None:
        bad_seqs = []
    return (pd.read_csv(f, sep='\s+', header=None)
     .rename(columns={0:'count', 1:'seq'})
     .query('count > @threshold')
     .query('seq != @bad_seqs') # remove any provided bad sequences
     .assign(fraction=lambda x: x['count']/x['count'].sum())
     .assign(file=f)
     .assign(ID=lambda x: x['file'].str.extract('(L\d+)'))
     .assign(seq=lambda x: x['seq'].pipe(fix_prefix_suffix))
     .assign(length=lambda x: x['seq'].str.len())
     .assign(rank=lambda x: np.arange(len(x)).astype(int)) 
    )

def make_wide_table(files):
    """ files = glob('Figure_2/TM_hists/*hist')
    """ 
    df = pd.concat([load_hist(f) for f in files])

    df_wide = df.pivot_table(index='seq', columns='ID', 
                             values='count')

    columns = natsorted(df_wide.columns)
    df_wide = df_wide[columns]
    df_wide['length'] = [len(s) for s in df_wide.index]
    return df_wide

def get_library(s):
    re.findall('L\d+', s)[0]

def fix_prefix_suffix(sequences, prefix='CACCG', suffix='GTTTT'):
    # checks the first string
    first_string = list(sequences)[0]
    m, n = len(prefix), len(suffix)
    if first_string.startswith(prefix) and first_string.endswith(suffix):
        # s[m:-n] starts at end of prefix and stops at beginning of suffix
        return [s[m:-n] for s in sequences]
    else:
        return sequences

def get_top_bcs(df_wide, n):
    top_bcs = set()
    for col in df_wide:
        top_bcs |= set(df_wide[col].sort_values(ascending=False)[:n].index)
    return top_bcs

def pairwise_distances(barcodes):
    """Returns pairwise distance matrix as a dataframe with labeled
    index (rows) and columns
    """
    n = len(barcodes)
    
    if n > 10000:
        print ('too many barcodes')
        return

    arr = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            arr[i,j] = arr[j, i] = distance(barcodes[i], barcodes[j])
            
    return pd.DataFrame(arr, index=barcodes, columns=barcodes)

def clustermap_from_distance_matrix(df_distances):
    """Create seaborn clustermap from custom distance matrix. Input the distance
    matrix as a dataframe to get labeled index (rows) and columns.
    """
    import scipy.spatial as sp
    import scipy.cluster.hierarchy as hc
    linkage = hc.linkage(sp.distance.squareform(df_distances), method='average')
    cg = sns.clustermap(df_distances, row_linkage=linkage, col_linkage=linkage)
    return cg

def add_rank(df_wide, rank):
    df_wide = df_wide.copy()
    drugs = ['JQ1', 'DMSO']
    for drug in drugs:
        sorter = lambda df: np.sort(df.values, axis=1)
        col = '%s_r%d' % (drug, rank)
#         col = (drug, rank)
        sorted_values = df_wide.fillna(-10).filter(regex='%s-\d' % drug).pipe(sorter)
#         sorted_values = df_wide.fillna(-10).pipe(sorter) 
#        print (sorted_values)
        df_wide[col] = sorted_values[:, rank]
#         print df_wide.fillna(-10).filter(regex='%s_\d' % drug).pipe(sorter).shape
        #print df_wide.fillna(-10).pipe(sorter).shape
    return df_wide

# bad_barcodes = pd.read_csv('bad_barcodes.csv', header=None)[0].pipe(list)
def load_TM_hist(f, bad_seqs=None, threshold=3):
    if bad_seqs is None:
        bad_seqs = []
    return (load_hist(f)
     .assign(ID=lambda x: x['file'].str.extract('(L\d+)'))
     .assign(seq=lambda x: x['seq'].pipe(fix_prefix_suffix))
#     .assign(log10=lambda x: x['fraction'].pipe(np.log10)) 
     .assign(log10=lambda x: np.log10(x['fraction']))
    .query('seq != @bad_seqs')
     .assign(length=lambda x: x['seq'].str.len())
)
def get_top_bcs(df_wide, n):
    top_bcs = set()
    for col in df_wide:
        top_bcs |= set(df_wide[col].sort_values(ascending=False)[:n].index)
    return top_bcs

bad_reps = ['d52, DMSO-4']  