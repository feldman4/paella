import functools
import pandas as pd


def bin_join(xs, symbol):
    symbol = ' ' + symbol + ' ' 
    return symbol.join('(%s)' % x for x in xs)
        
or_join  = functools.partial(bin_join, symbol='|')
and_join = functools.partial(bin_join, symbol='&')


def set_apply(df, col, f):
	"""Join the results of applying `f` (returns a dictionary).
	"""
	inputs = sorted(set(df[col]))
	result = pd.DataFrame([f(x) for x in inputs], index=inputs)
	return df.join(result, on=col)

def transform_cols(df, regex, f):
    df = df.copy()
    cols = df.filter(regex=regex).columns
    df[cols] = f(df[cols])
    return df


def add_statistic(df, col_regex, col_name, statistic):
    """Calculate a statistic some columns. 
    """
    df_ = df.filter(regex=col_regex)
    if statistic == 'mean':
        stat = df_.mean(axis=1).rename(col_name + '_' + statistic)
    elif statistic == 'min':
        stat = df_.min(axis=1).rename(col_name + '_' + statistic)
    elif statistic == 'max':
        stat = df_.max(axis=1).rename(col_name + '_' + statistic)
    else:
        stat = statistic(df_).rename(col_name + '_' + 'custom')
        
    return pd.concat([df, stat], axis=1)