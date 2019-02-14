import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import paella.utils
import seaborn as sns

D458_examples = ['L03', 'L04', 'L05', 'L26', 'L27']

D458_ETP = ['L26'] 
D458_JQ1 = ['L27', 'L28', 'L29', 'L30', 'L31']
D458_DMSO = ['L32', 'L33', 'L34', 'L35', 'L36']
D458_DMSO_1 = ['L32', 'L33', 'L34', 'L36'] # drop bad replicate

D458_screen = D458_ETP + D458_JQ1 + D458_DMSO

PC9_ETP = ['L10']
PC9_1uM  = ['L11', 'L12', 'L13', 'L14', 'L15']
PC9_60nM = ['L16', 'L17', 'L18', 'L19', 'L20']
PC9_DMSO = ['L21', 'L22', 'L23', 'L24', 'L25']

PC9_screen = PC9_ETP  + PC9_60nM + PC9_1uM + PC9_DMSO
PC9_screen_1uM = PC9_ETP  + PC9_1uM + PC9_DMSO

PC9_ETP_20N = ['L178', 'L179', 'L180', 'L181', 'L182']
PC9_ETP_20N_1 = ['L179']
PC9_1uM_20N = ['L169', 'L170', 'L171', 'L172', 'L173']
PC9_Osi = ['L174', 'L175', 'L176', 'L177']
PC9_DMSO_20N = ['L164', 'L165', 'L166', 'L167', 'L168']

PC9_ETP_20N = ['L178', 'L179', 'L180', 'L181', 'L182']


def plot_count_dist(counts, ticks=None):

    fig, ax0 = plt.subplots()
    ax1 = ax0.twinx()

    bins = list(counts.index)

    if ticks is None:
        ticks = [min(bins)] + [10, 30, 100, 300, 1000]


    cdf = np.cumsum(counts) / sum(counts)
    ax1.plot(np.log10(bins), cdf)
    ax1.set_ylabel('cumulative fraction')
    ax0.set_ylabel('number of barcodes')

    ax0.stem(np.log10(bins), counts, markerfmt='.')

    
    ax0.set_xticks([np.log10(x) for x in ticks])
    ax0.set_xticklabels(ticks);

    # plt.tight_layout(h_pad=-100)

    return fig, ax0, ax1


def plot_D458_abundances(df_wide, df_info):
    fig, ax = plt.subplots()

    for library in D458_examples:
        xs = df_wide[library].dropna().sort_values(ascending=False)
        label = df_info.loc[library]['name']
        ax.plot(xs[1:200000].pipe(list), label=label)
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([1e2, 1e6])    
    ax.legend(bbox_to_anchor=(1.07, 1.03))
    return fig


def load_PC9_cell_counts(f):
    columns = {0: 'count', 'column': 'mimi_did', 'Day': 'day', 'level_0': 'sample'}
    pat = '(.*?)_?\d$'
    
    df_counts = (pd.read_csv(f, index_col=0, skiprows=1, header=[0,1]).iloc[:16, :12]
     .stack().stack().reset_index()
     .rename(columns=columns)
     .assign(day=lambda x: x['day'].astype(int))
     .assign(drug=lambda x: x['sample'].str.extract(pat)[0])
    )
    
    df_counts_wide = (df_counts
     .query('mimi_did != "frozen"')
     .pivot_table(index=['sample', 'drug'], columns=['mimi_did', 'day'], values='count')
    )

    return df_counts, df_counts_wide


def PC9_growth_rate(df_counts_wide):
    # make a growth rate table (doublings per day)
    growth_rate = (df_counts_wide['Cells counted'].values / 
                   df_counts_wide['Cells plated'].values)

    missing_data = (growth_rate == 0 | np.isnan(growth_rate))
    growth_rate[missing_data] = np.nan

    intervals = list(df_counts_wide['Cells counted'].columns - 
                     df_counts_wide['Cells plated'].columns)

    data = growth_rate / intervals

    x = df_counts_wide['Cells counted']
    df_growth_wide = pd.DataFrame(data, index=x.index, columns=x.columns)
    df_growth = (df_growth_wide
     .stack().reset_index().rename(columns={0: 'growth_rate'})
     )

    return df_growth, df_growth_wide


def PC9_cumulative_count(df_counts, df_counts_wide, ):

    # use plating ratio to adjust cell counts
    plating_days = [3, 7, 12, 17]
    plating_ratio = (df_counts_wide['Cells counted'][plating_days] / 
                     df_counts_wide['Cells plated'][plating_days])


    df_adj_wide = df_counts_wide['Cells counted'].copy()

    offset_days = [7, 12, 17, 22]
    df_adj_wide[offset_days] *= plating_ratio.cumprod(axis=1).fillna(1).values

    df_adj = (df_adj_wide
     .stack().reset_index().rename(columns={0: 'adj_count'}))

    df_day1 = (df_counts.query('day == 1')
               .rename(columns={'count': 'adj_count'}))

    df_adj = pd.concat([df_day1,
                   df_adj], sort=True)

    return df_adj


def plot_barcode_color(df_plt, hue_order, palette, optimal_ordering=True):
    """Plot clustermap using `barcode_set` column to label rows.
    """
    assert df_plt['barcode_set'].pipe(set) == set(hue_order)
    get_color = lambda x: dict(zip(hue_order, palette)).get(x, 'white')
    colors = df_plt['barcode_set'].map(get_color)
    
    df_plt = df_plt.drop(['barcode_set'], axis=1) 

    from scipy.cluster import hierarchy
    row_linkage = hierarchy.linkage(df_plt, optimal_ordering=optimal_ordering)

    return sns.clustermap(data=df_plt, 
        #row_linkage=row_linkage,
        row_cluster=False,
        row_colors=colors,  col_cluster=False)
     
    

def plot_stacked_replicates(df_wide, threshold, colors=None, reverse_legend=True):
    df = (df_wide
      .pipe(lambda x: x>threshold)
      .assign(count=lambda x: x.sum(axis=1))
      )
    
    arr = []
    for col in df_wide.columns:
        # only take barcodes present in this JQ1 sample
        filt = df[col] > 0
        # for those barcodes, histogram of number of JQ1 replicates with the barcode
        counts = df[filt]['count'].value_counts().rename(col)
        counts.index + 1
        arr.append(counts)
    
    df_shared_reps = pd.concat(arr, axis=1)
    
    if colors is None:
        colors = sns.color_palette('Reds', df_wide.shape[1] - 1) + ['black']
    ax = (df_shared_reps
     .pipe(lambda x: x / x.sum())
     [::-1] #reverse d
     .T.plot(kind='bar', stacked=True, 
             color=colors,
              grid = False,
            )
    )

    ax.set_ylabel('fraction of barcodes in sample')
    ax.set_xticklabels(ax.get_xticklabels(), rotation='horizontal')
    
    handles, labels = ax.get_legend_handles_labels()
    if reverse_legend:
        handles, labels = handles[::-1], labels[::-1]
    ax.legend(handles, labels, bbox_to_anchor=(1, 0.85), frameon=False)
        
    return df_shared_reps, ax
