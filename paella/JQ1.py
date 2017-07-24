from paella.imports import *

from scipy.stats import pearsonr
pearsonr_ = lambda a, b: pearsonr(a, b)[0]

def load_reads():
    f0 = 'JQ1_60K/207KBarcodes_d458.gct.txt'
    f1 = 'JQ1_60K/scores-Mimi_600KBC_20150113_D283.gct.txt'
    reads0 = (pd.read_csv(f0, sep='\t', skiprows=2)
                .set_index('Name')
                .rename(columns={'Description': 'Sum'}))

    reads0.columns = ['D458_' + c for c in reads0.columns]

    reads1 = (pd.read_csv(f1, sep='\t', skiprows=2)
                .rename(columns={'Construct Barcode': 'Name'})
                .set_index('Name')
                .drop('Construct IDs', axis=1))

    reads = pd.concat([reads0, reads1], axis=1).fillna(0).astype(int)
    
    return reads

def filter_D458(reads, thresholds=(3,1000)):
    low, high = thresholds
    exclude = ['D458_JQ1_5']
    filt = (reads['D458_ETP'] > low) & (reads['D458_ETP'] < high)
    return reads[filt].drop(exclude, axis=1).filter(like='D458')

def filter_D283(reads, thresholds=(3,1000)):
    low, high = thresholds
    exclude = [c for c in reads.columns if 'FT' in c or 'JQ1_5' in c]
    filt = (reads['D283_ETP'] > low) & (reads['D283_ETP'] < high) 
    return reads[filt].drop(exclude, axis=1).filter(like='D283')


def add_shared(reads_, col_filter, threshold=3):
    reads_ = reads_.copy()
    cols = [c for c in reads_.columns if col_filter in c]
    shared = 'shared_%s' % col_filter
    reads_[shared] = (reads_ > threshold)[cols].sum(axis=1)
    return reads_

def plot_shared(reads_, col_filter, threshold=3):
    shared = 'shared_%s' % col_filter
    filt = reads_[shared] > 0
    
    bins = np.arange(10) - 0.5
    ax = (reads_.loc[filt, shared] - 1).hist(bins=bins, width=0.7)
    ax.set_title('# shared in %s, reads > %d' % (col_filter, threshold))
    return ax

def plot_stacked_bar_JQ1(reads_, ax):
    colors = 'black', '#fee5d9','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#99000d'

    arr = []
    JQ1_cols = [c for c in reads_.columns if 'JQ1_' in c]
    for i, col in enumerate(JQ1_cols):
        filt = reads_[col] > 3
        cumulative = (reads_.loc[filt, 'shared_JQ1']
                            .value_counts()
                            .sort_index(ascending=False)
                            .cumsum())
        cumulative = 100. * (cumulative / float(cumulative.max()))
        for j, (y, color) in enumerate(zip(cumulative[::-1], colors)):
            label = str(j) if i == 0 else None
            ax.bar(i, y, color=color, label=label)


    ax.set_xticks(range(len(JQ1_cols)))
    ax.set_xticklabels(JQ1_cols, rotation=90)

    legend = ax.legend(title='# shared', bbox_to_anchor=(1., 1.))
    plt.setp(legend.get_title(),fontsize=14)
    return


def plot_correlation(reads_, cell_line):
#     cell_line = 'D458'
    drop = reads_.filter(regex='Sum|shared').columns
    corr = pairwise.compare(reads_.drop(drop, axis=1), pearsonr_)
    
    fig, ax = plt.subplots()
    ax = sns.heatmap(corr)
    ax.set_title('pearson correlation')
    ax.set_xlabel(''); ax.set_ylabel('')
    
    cg = sns.clustermap(corr)
    cg.ax_heatmap.set_title('%s pearson correlation' % cell_line)
    cg.ax_col_dendrogram.set_visible(False)
    cg.ax_heatmap.yaxis.set_visible(False)
    cg.ax_heatmap.set_xlabel('')

    return fig, cg, corr