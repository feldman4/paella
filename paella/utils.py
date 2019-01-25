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