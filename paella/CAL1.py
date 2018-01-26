from paella.imports import *
from paella.data import *
from glob import glob

fig_path = 'figures/20170816_CAL-1/'
files_0 = glob('NGS/20170816/hist/*hist')
plasmid_file = '/Users/feldman/paella/20170213_DF/hist/C04.grep.hist'
files = list(files_0) + [plasmid_file]

def load():
	pre_filter = 3
	arr = []
	for f in files:
	    if 'C04' in f:
	        sample = 'pDNA'
	    else:
	        sample = re.findall('(L..)', f)[0]
	    
	    df = read_hist(f, strip=5)
	    if pre_filter:
	        df = df[df['count'] > pre_filter]
	    df['sample'] = samples[sample]
	    df['rank'] = df.index + 1
	    arr += [df]

	df = pd.concat(arr)
	df = df.reset_index(drop=True)
	df['reads_in_sample'] = df.groupby('sample')['count'].transform(sum)
	df['fraction'] = df['count'] / df['reads_in_sample']


	filt = df['count'] > 3
	mapped_reads = df.groupby('sample')['reads_in_sample'].apply(max)
	barcodes_gt3 = df[filt].groupby('sample').apply(len)
	df_info = pd.DataFrame([mapped_reads, barcodes_gt3]).T
	df_info.columns = 'mapped reads', 'barcodes with 3+ reads'
	df_info
	return df, df_info


def add_resampling(df):