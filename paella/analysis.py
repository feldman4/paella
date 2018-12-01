from itertools import product
from paella.imports import *

def rc(seq):
    return ''.join(watson_crick[x] for x in seq)[::-1]

def degenerate(s):
    """a generator would be nice
    """
    bases = {'N': list('ACTG')}
    arr = []
    for c in s:
        arr += [bases.get(c, [c])]
    return [''.join(x) for x in product(*arr)]


watson_crick = {'A': 'T',
                'T': 'A',
                'C': 'G',
                'G': 'C',
                'U': 'A',
                'N': 'N'}
watson_crick.update({k.lower(): v.lower() for k, v in watson_crick.items()})

def read_fastq(f):
    with open(f, 'r') as fh:
        return fh.read().split('\n')[1::4]

def read_fasta(f):
    """fasta file with one sequence
    """
    with open(f, 'r') as fh:
        lines = fh.read().split('\n')
    header = lines[0][1:]
    seq = ''.join(lines[1:])
    return header, seq

def to_96W(df, values):
    return df.pivot(index='row', columns='col', values=values).fillna(0)

def file_table(fastq_glob):
    """ 
    Read counts from fastq files. 
    Ex: file_table('fastq/[A-H]*R1*.fastq')
    """

    # file_info = !wc -l fastq/[A-H]*R1*.fastq
    file_info = IPython.get_ipython().system('wc -l ' + fastq_glob)
    df_files = pd.DataFrame([f.split() for f in file_info[:-1]]
                            , columns=['reads', 'fastq'])
    df_files['reads'] = [int(x)/4   for x in df_files['reads']]
    df_files['well']  = [s.split('/')[1][:3] for s in df_files['fastq']]
    df_files['row']   = [w[0]       for w in df_files['well']]
    df_files['col']   = [int(w[1:]) for w in df_files['well']]

    return df_files

def read_hist(f, strip=0):
    """
    Read output of sort | hist -c | sort -bnr.
    """
    df = pd.read_csv(f, sep='\s+', header=None)
    df.columns='count', 'seq'
    if strip:
        df['seq'] = [s[strip:-strip] for s in df['seq']]

    return df

def make_well_pairs(wells):
    x = list(combinations(wells, 2))
    return pd.DataFrame(x, columns = ['A', 'B'])


def calculate_diversity(df, n):

    arr = [list(s) for s in df['sgRNA'][:n]]
    dfc = pd.DataFrame(arr)
    dfc = dfc.stack().reset_index()
    dfc.columns = 'barcode', 'position', 'base'
    dfc = dfc[dfc['base'] != 'N']

    from scipy.stats import entropy
    qits = []
    for i in range(26):
        df = dfc[dfc['position']==i]
        p = Counter(df['base']).values()
        qits += [4**entropy(p, base=4)] 

    return dfc, qits

def stop_codon(sequence, n=20):
    """Provide an sgRNA ending at (not including) GTTT.
    """
    sequence = sequence[-n:].upper()
    return 'ATG' in sequence and 'CAT' in sequence 


# def apply_overlap(df_overlap, df_counts, f):
#   for wellA, wellB in zip(df_overlap[])


def resample(df_, n=350000):
    x = df_['count'].copy()
    x /= x.sum()
    c = Counter(np.random.choice(x.index, size=n, p=x))
    s = pd.Series(c.values(), index=c.keys())
    s.name = 'resampled_count'
    df2 = pd.concat([x, s], axis=1).fillna(0).astype(int)
    return df2['resampled_count']

def well_to_row_col(df, in_place=False, col_to_int=True):
    if not in_place:
        df = df.copy()
    f = lambda s: int(s) if col_to_int else s

    df['row'] = [s[0] for s in df['well']]
    df['col'] = [f(s[1:]) for s in df['well']]

    return df