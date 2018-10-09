from paella.imports import *

def strip_riz(s):
    """68 bp P5 read not enough
    """
    front = 'AAACACCG'
    back  = 'GT'
    return s[len(front):-1*len(back)]

def strip_sg(s):
    front = 'ACCG'
    back = 'GTTT'
    return s[len(front):-1*len(back)]

def load_hist(f):
    df = pd.read_csv(f, sep='\s+', header=None)
    df.columns = 'count', 'seq'
    df['well'] = get_sample(f)
    return df

def get_sample(f):
    pat = 'G_(...)_'
    return re.findall(pat, f)[0]

def calc_fraction(df):
    return df['count'] / df['count'].sum()

def calc_rank(df):
    # first argsort returns indices that would sort the list
    return np.argsort(np.argsort(-1 * df['count']))

def fraction_and_rank(df):
    return (# apply `calc_fraction`, name the resulting column "fraction"
            df.assign(fraction=calc_fraction)
              .assign(rank=calc_rank))

def min_rank(x):
    return x.filter(regex='rank', axis=1).min(axis=1)    

def pivot_fraction_rank(df, fill_rank=None, fill_frac=None):
    """Pivot based on `seq`, concatenating fraction and rank. Optionally fill in 
    na values. If `fill_rank` is given, the ranks are returned as `int` type.
    """
    df_rank = df.pivot_table(values='rank', index='seq', columns='sample_info')
    if fill_rank is not None:
         df_rank = (df_rank.fillna(fill_rank)
                          .astype(int))

    # maintain original sample order after pivot
    sample_order = list(df['sample_info'].drop_duplicates())
    df2 = (df.pivot_table(values='fraction', index='seq', columns='sample_info')
    	     [sample_order]
             .join(df_rank, rsuffix='_rank', lsuffix='_frac')
             .assign(min_rank=min_rank))
    
    if fill_frac is not None:
        return df2.fillna(1e-6)
    return df2
