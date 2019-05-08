import numpy as np
import pandas as pd

import paella.utils
import seaborn as sns

D458_examples = ['L32', 'L26', 'L27','L05']
D458_examples = [
    ('L26', 'ETP'), 
    ('L32', 'DMSO'),
    ('L27', 'JQ1')]  

D458_all_examples = {'DMSO-1': 'L32',
 'DMSO-2': 'L33',
 'DMSO-3': 'L34',
 # 'DMSO-4': 'L35',
 'DMSO-5': 'L36',
 'ETP': 'L26',
 # 'JQ1-1': 'L27',
 'JQ1-2': 'L28',
 'JQ1-3': 'L29',
 'JQ1-4': 'L30',
 'JQ1-5': 'L31'}

#D458 26N
D458_pre_freeza = ['L03']
D458_post_freeze= ['L04']
D458_ETP = ['L26'] 
D458_JQ1 = ['L27', 'L28', 'L29', 'L30', 'L31']
D458_DMSO = ['L32', 'L33', 'L34', 'L35', 'L36']
D458_DMSO_1 = ['L32', 'L33', 'L34', 'L36'] # drop bad replicate

D458_screen = D458_ETP + D458_JQ1 + D458_DMSO

#PC9 26N
PC9_ETP = ['L10']
PC9_1uM  = ['L11', 'L12', 'L13', 'L14', 'L15']
PC9_60nM = ['L16', 'L17', 'L18', 'L19', 'L20']
PC9_DMSO = ['L21', 'L22', 'L23', 'L24', 'L25']

PC9_screen = PC9_ETP  + PC9_60nM + PC9_1uM + PC9_DMSO
PC9_screen_1uM = PC9_ETP  + PC9_1uM + PC9_DMSO

#PC9 20N
PC9_ETP_20N = ['L178', 'L179', 'L180', 'L181', 'L182']
PC9_ETP_20N_1 = ['L179']
PC9_1uM_20N = ['L169', 'L170', 'L171', 'L172', 'L173']
PC9_Osi = ['L174', 'L175', 'L176', 'L177']
PC9_DMSO_20N = ['L164', 'L165', 'L166', 'L167', 'L168']

#Hela 20N
Hela_ETP = ['L218']
Hela_PBS = ['L219', 'L220', 'L221', 'L222', 'L223']
Hela_Hygro = ['L224', 'L225', 'L226', 'L227', 'L228']

T1 = 'AGATCGTACCAGGGATTGGG'
T2 = 'TGCAGTGGCCGTTGACAAAT'
T3 = 'GCGGGATCATTGCAATTATA'


custom_rcParams = {
    'figure.figsize':(9, 8),
    'legend.frameon': False,
    'axes.labelsize' : 40,
    'xtick.labelsize': 40,
    'ytick.labelsize': 40,
}

def apply_rcParams():
    from matplotlib import rcParams
    for k in custom_rcParams:
        rcParams[k] = custom_rcParams[k]


def plot_count_dist(counts, show_total=True, ticks=None):

    fig, ax0 = plt.subplots()
    ax1 = ax0.twinx()

    bins = list(counts.index)

    if ticks is None:
        ticks = [min(bins)] + [10, 30, 100, 300, 1000]


    cdf = np.cumsum(counts) / sum(counts)
    ax1.plot(np.log10(bins), cdf)
    ax1.set_ylabel('Cumulative fraction')
    ax0.set_ylabel('Number of barcodes')
    ax0.set_xlabel('Read counts')

    ax0.stem(np.log10(bins), counts, markerfmt='.')

    
    ax0.set_xticks([np.log10(x) for x in ticks])
    ax0.set_xticklabels(ticks);

    # plt.tight_layout(h_pad=-100)

    if show_total:
        label = '{:,} barcodes detected'.format(int(counts.sum()))
        ax0.text(0.32, 0.75, label, transform=ax0.transAxes, fontsize=16)
    return fig, ax0, ax1


def plot_D458_abundances(df_wide, df_info):
    fig, ax = plt.subplots()

    for library in D458_examples:
        xs = df_wide[library].dropna().sort_values(ascending=False)
        # label = df_info.loc[library]['name']
        label = df_info.loc[library]['sample_info']
        ax.plot(xs[1:200000].pipe(list), label=label)
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    # ax.set_xlim([1e2, 1e6]) 
    ax.set_xlim([1e0, 1e6])    
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
        cmap='RdBu_r',
        row_cluster=False,
        row_colors=colors,  
        yticklabels=0,
        xticklabels=0,
        col_cluster=False)
     
    

def plot_stacked_replicates(df_wide, threshold, colors=None, reverse_legend=False):
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

    ax.set_ylabel('Fraction of JQ1-enriched \nbarcodes in replicates')
    ax.set_xlabel('JQ1 replicates')
    ax.set_xticklabels(ax.get_xticklabels(), rotation='horizontal')
    # legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #    ncol=2, mode="expand", borderaxespad=0.)
    
    handles, labels = ax.get_legend_handles_labels()
    if reverse_legend:
        handles, labels = handles[::-1], labels[::-1]
    ax.legend(handles, labels, bbox_to_anchor=(0.0, 1.02, 1.1, 0.102), 
        frameon=False, 
        loc=6, ncol=5, mode="expand", borderaxespad=0.1, fontsize=30)
    # ax.legend(handles, labels, bbox_to_anchor=(-0.2, 1.02, 1.5, 0.102), 
    #     frameon=False, 
    #     loc=6, ncol=5, mode="expand", borderaxespad=0.009, fontsize=20)
        
    return df_shared_reps, ax

def log_scale_ticks(ax, axis='x'):
    """Format x tick labels as exponents
    """
    if axis == 'x':
        get_labels, set_labels = ax.get_xticks, ax.set_xticklabels
    elif axis == 'y':
        get_labels, set_labels = ax.get_yticks, ax.set_yticklabels
    else:
        raise ValueError
    labels = get_labels()
    format_exponent = lambda x: '{' + '{0:d}'.format(int(x)) + '}'
    new_labels = ['$10^{0}$'.format(format_exponent(x)) for x in labels]
    set_labels(new_labels)
