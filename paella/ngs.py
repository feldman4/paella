import numpy as np
import pandas as pd
from natsort import natsorted
from Levenshtein import distance
import re

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
     .assign(ID=re.findall('L\d+', f)[0])
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



def pairwise_distances(barcodes):
    """Returns pairwise distance matrix as a dataframe with labeled
    index (rows) and columns
    """
    n = len(barcodes)
    
    if n > 10000:
        print('too many barcodes')
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


def cleanup_duplicates(df, top_n, convert_log=True, kmer=None):
    from paella.UMIGraph import UMIGraph
    """Sort input, duplicates only removed from top n.
    Discards duplicate entries (could also sum).
    """ 
    cols = list(df.filter(regex='^L').columns)
    data = df[cols]
    if convert_log:
        data = 10**(7 + data)

    # the actual count used for each barcode is just averaged across columns
    # in the table
    counts = data.mean(axis=1)

    G = UMIGraph(df.index[:top_n], 
                 counts=counts[:top_n].astype(int), 
                 threshold=3, kmer=kmer)
    discard = sorted(set(df.index[:top_n]) - set(G.component_dict.values()))
    print('discarding', len(discard), 'duplicates')
    return df.query('index != @discard')

def fillna_wide(df, log_cutoff):
    """Use different rules for libraries (starts with L) and ranks.
    """
    df = df.copy()
    cols  = [x for x in df.columns if x.startswith('L')]
    cols_ = [x for x in df.columns if 'rank' in x]
    df[cols] = df[cols].fillna(log_cutoff)
    df[cols_] = df[cols_].fillna(df[cols_].max().max())
    return df


def get_rank(df_wide, rank):
    return (df_wide
     .pipe(lambda x: (-x).rank())
     .pipe(lambda x: np.sort(x)[:, rank])
    )

def get_val(df_wide, rank):
    """rank is high to low (e.g., rank=0 returns max)
    """
    return -np.sort(-df_wide.values, axis=1)[:, rank]


def get_sample_groups(df_wide, groups, log_cutoff, ranks=('max', '2nd', 'med')):
    """Get subset of wide table with samples in groups, assign ranks within
    groups. 
    groups = {'JQ1': figure2.D458_JQ1, 'DMSO': figure2.D458_DMSO}
    """
    def add_ranks(df, name, ranks):
        d = {}
        for r in ranks:
            rank_name = 'rank_{0}_{1}'.format(r, name)
            val_name = 'val_{0}_{1}'.format(r, name)
            if r == 'max':
                d[rank_name] = get_rank(df, 0)
                d[val_name]  = get_val(df, 0)
            elif r == '2nd':
                # 2nd lowest
                if df.shape[1] > 1:
                    d[rank_name] = get_rank(df, -2)
                    d[val_name]  = get_val(df, -2)
            elif r == 'med':
                n = df.shape[1]
                if 2 * int(n/2) == n:
                    # even
                    i0, i1 = int(n/2) - 1, int(n/2)
                    d[rank_name] = (get_rank(df, i0) + get_rank(df, i1))/2
                    d[val_name]  = (get_val(df, i0)  + get_val(df, i1))/2
                else:
                    # odd
                    i =  int(n/2)
                    d[rank_name] = get_rank(df, i)
                    d[val_name] = get_val(df, i)
            else:
                raise ValueError('{0} not recognized'.format(r))
        return df.assign(**d)
    arr = []
    for name, samples in groups.items():
        (df_wide
         [samples]
         .dropna(how='all').fillna(log_cutoff)
         .pipe(add_ranks, name, ranks)
         .pipe(arr.append)
        )
    return pd.concat(arr, axis=1, sort=True).pipe(fillna_wide, log_cutoff)


def assign_barcode_sets(df_wide, group_names, num_top, prefix='rank_max_'):
    """Define barcode sets based on rank threshold in each group (e.g.,
    set of top barcodes in one group only, set of top barcodes across all 
    groups, etc). Applies threshold to columns with name '{prefix}{group}' (e.g., 
    'rank_min_DMSO').

    `num_top` can be a single number or a dictionary from group names to thresholds.
    """
    # assign barcode sets
    cols = [prefix + c for c in group_names]
    if not isinstance(num_top, dict):
        num_top = {k: num_top for k in group_names}

    arr = []
    for vals in df_wide[cols].values:
        xs = []
        for c, v in zip(group_names, vals):
            if v < num_top[c]:
                xs += ['{0}_{1}'.format(c, num_top[c])]
        arr.append('_'.join(xs))

    return df_wide.assign(barcode_set=arr)