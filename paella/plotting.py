from paella.imports import *


def plot_psn(lrange=(0.1, 0.3, 1), krange=tuple(range(1,10))):
    from math import factorial
    psn = lambda l, k: (l**k / factorial(k)) * np.exp(-l)
    lrange = [0.1, 0.3, 1]
    krange = range(1,10)
    arr = []
    for l in lrange:
        for k in krange:
            arr += [[l, k, psn(l, k)]]    

    df_ = pd.DataFrame(arr, columns=('l', 'k', 'count'))
    df_ = df_.pivot_table(index='k', columns='l', values='count')
    df_ = df_ / df_.sum()
    
    return df_.plot()

def plot_reg_edit_distance(df):
    df = df.copy()
    df['1'] = df.loc[:,1]
    ax = sns.regplot(data=df, x='count', y='1', x_estimator=np.mean)
    ax.set_ylabel('sequencing errors\n(barcodes with edit distance=1)')
    ax.set_xlabel('barcode read count')
    ax.set_title('reads attributable to sequencing error vs. barcode abundance')
    ax.figure.tight_layout()
    return ax

def add_distance_counts(df, k=26, n=100):
    df = df.copy()
    for j in range(k):
        df[j] = np.nan

    n = 1000
    arr = []
    all_distances = []
    for i in range(n):
        s = df.ix[i, 'sgRNA']
        distances = [distance(s, t) for t in df['sgRNA']]
        arr += [np.bincount(distances, minlength=k)]
        all_distances += [distances]

    df.loc[:n-1,range(k)] = np.array([x[:k] for x in arr])
    
    return df, all_distances

def plot_cdf(counts, ax):
    cs  = np.cumsum(counts)
    cs /= float(sum(counts))
    xr  = np.arange(len(counts)) / float(len(counts))
    ax.plot(xr, cs)
    return xr, cs

def plot_hist(df, norm=True, axs=None, range_=None, zoom_xlim=10, cdf_cutoff=3):
    if axs is None:
        _, (ax0, ax1, ax2) = plt.subplots(1,3)
    else:
        ax0, ax1, ax2 = axs
        
    if range_ is None:
        range_ = range(df['count'].max())
       
    vals = np.bincount(df['count'])
    ax0.plot(vals); ax0.set_yscale('log'); ax0.set_xscale('log')
    ax0.set_xlim([0, 1e4])
    df['count'].hist(bins=range(11), ax=ax1).set_yscale('log')
    ax1.set_xlim([0, zoom_xlim]); ax1.get_yaxis().set_visible(False)

    cdf_counts = df['count']
    cdf_counts = cdf_counts[cdf_counts > cdf_cutoff] # read cutoff
    plot_cdf(cdf_counts, ax2)
    ax2.get_yaxis().set_visible(False)

    ax0.set_ylim([1, 1e7])
    ax1.set_ylim([1, 1e7])
    return ax0, ax1, ax2


### EXPERIMENTS

def plot_rank_abundance_hela(df):
    # HeLa
    df2 = df.query('rank < 150')

    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_yscale('log')
    style = { '1e4_Hygro':      dict(color='blue', linestyle='-')
            , '3.5e5_Hygro':    dict(color='blue', linestyle=':')
            , '1e4_PBS':        dict(color='gray', linestyle='-')
            , '3.5e5_PBS':      dict(color='gray', linestyle=':')
            , 'HeLa_1e4_ETP':   dict(color='green', linestyle='-')
            , 'HeLa_3.5e5_ETP': dict(color='green', linestyle=':')}
    for sample, df_ in df2.groupby('sample'):
        try:
            ax.plot(df_['fraction'], label=sample[:-1], **style[sample[5:-1]])
        except KeyError:
            ax.plot(df_['fraction'], label=sample, **style[sample])
            
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())
    ax.set_xlabel('rank')
    ax.set_ylabel('abundance')
    return ax

def plot_correlation_hela(df, bottleneck):

    color_top, color_mid, color_mid_shared = 10, 10, 10
    threshold, mid_shared_threshold = 2e-5, 4
    mid = 100
    ds = 1

    bottlenecks = {
        '4e4':
                   ['HeLa_4e4_ETP'
                   ,'HeLa_4e4_PBS1'
                   ,'HeLa_4e4_PBS2'
                   ,'HeLa_4e4_PBS3'
                   ,'HeLa_4e4_PBS4'
                   ,'HeLa_4e4_Hygro1'
                   ,'HeLa_4e4_Hygro2'
                   ,'HeLa_4e4_Hygro3'
                   ,'HeLa_4e4_Hygro4'
                   ]
        ,'3.5e5': 
               ['HeLa_3.5e5_ETP'
               ,'HeLa_3.5e5_PBS1'
               ,'HeLa_3.5e5_PBS2'
               ,'HeLa_3.5e5_PBS3'
               ,'HeLa_3.5e5_PBS4'
               ,'HeLa_3.5e5_Hygro1'
               ,'HeLa_3.5e5_Hygro2'
               ,'HeLa_3.5e5_Hygro3'
               ,'HeLa_3.5e5_Hygro4'
               ]
           }

    samples = bottlenecks[bottleneck]

    filt  = df['sample'].isin(samples)

    df_plot = (df[filt].pivot(index='sgRNA', columns='sample', values='fraction')
                   .iloc[::ds]).copy()

    keep_filt = ~((df_plot.isnull()).sum(axis=1) >= mid_shared_threshold)
    keep_filt = keep_filt[keep_filt].index
    df_plot[df_plot < threshold] = threshold
    df_plot = df_plot.fillna(threshold / 2)

    # color rank subgroups
    ranked = df_plot.filter(regex='Hygro').max(axis=1).sort_values(ascending=False)
    top = ranked[:color_top].index
    mid = ranked[mid:mid+color_mid].index
    mid_shared = ranked[ranked.index.isin(keep_filt)][100:100+color_mid_shared].index
    categories = 'log10 > %.1f' % np.log10(threshold), 'top', 'mid', 'mid_shared'
    palette = (0.5, 0.5, 0.5, 0.5), 'red', 'blue', 'green'
    df_plot['category'] = categories[0]
    df_plot.loc[top, 'category'] = categories[1]
    df_plot.loc[mid, 'category'] = categories[2]
    df_plot.loc[mid_shared, 'category'] = categories[3]

    # more filtering...
    highest = df_plot.filter(regex='PC9').max(axis=1)
    filt = (df_plot['category'] == categories[0]) & (highest <= threshold)
    print len(filt), filt.sum()

    float_cols = df_plot.filter(regex='HeLa').columns
    df_plot_ = df_plot[~filt].copy()

    A = (df_plot_[float_cols] <= threshold).as_matrix()
    B = (df_plot_['category'] == categories[0]).as_matrix()[:, None]

    data = df_plot_[float_cols].as_matrix()
    data[(A&B)] = np.nan

    df_plot_[float_cols] = data
    df_plot_[(df_plot_['category'] == categories[0])]
    df_plot_[float_cols] = np.log10(df_plot_[float_cols])
    pg = sns.pairplot(df_plot_, hue='category', palette=palette, hue_order=categories
                      , x_vars=samples, y_vars=samples, plot_kws={'s': 16})

    lim = [-5.2, 0]
    [(ax.set_xlim(lim), ax.set_ylim(lim)) for ax in pg.axes.flat[:]];
    # pg.fig.savefig('test.pdf')

    return ax

def plot_correlation_cal1(df, hela_too=False):

    color_top, color_mid, color_mid_shared = 10, 10, 10
    threshold, mid_shared_threshold = 2e-5, 4
    mid = 100
    ds = 1

    import paella.data
    filt  = df['sample'].apply(lambda x: 'Cal-1' in x)

    if hela_too:
        filt |= df['sample'].apply(lambda x: 'HeLa' in x)

    samples = sorted(set(df.loc[filt, 'sample']))

    df_plot = (df[filt].pivot(index='sgRNA', columns='sample', values='fraction')
                   .iloc[::ds]).copy()

    keep_filt = ~((df_plot.isnull()).sum(axis=1) >= mid_shared_threshold)
    keep_filt = keep_filt[keep_filt].index
    df_plot[df_plot < threshold] = threshold
    df_plot = df_plot.fillna(threshold / 2)

    # color rank subgroups
    ranked = df_plot.filter(regex='SL-401').max(axis=1).sort_values(ascending=False)
    top = ranked[:color_top].index
    mid = ranked[mid:mid+color_mid].index
    mid_shared = ranked[ranked.index.isin(keep_filt)][100:100+color_mid_shared].index
    categories = 'log10 > %.1f' % np.log10(threshold), 'top', 'mid', 'mid_shared'
    palette = (0.5, 0.5, 0.5, 0.5), 'red', 'blue', 'green'
    df_plot['category'] = categories[0]
    df_plot.loc[top, 'category'] = categories[1]
    df_plot.loc[mid, 'category'] = categories[2]
    df_plot.loc[mid_shared, 'category'] = categories[3]


    float_cols = df_plot.filter(regex='Cal-1').columns
    df_plot_ = df_plot.copy()

    A = (df_plot_[float_cols] <= threshold).as_matrix()
    B = (df_plot_['category'] == categories[0]).as_matrix()[:, None]

    data = df_plot_[float_cols].as_matrix()
    data[(A&B)] = np.nan

    df_plot_[float_cols] = data
    df_plot_[(df_plot_['category'] == categories[0])]
    df_plot_[float_cols] = np.log10(df_plot_[float_cols])
    pg = sns.pairplot(df_plot_, hue='category', palette=palette, hue_order=categories
                      , x_vars=samples, y_vars=samples, plot_kws={'s': 16})

    lim = [-5.2, 0]
    [(ax.set_xlim(lim), ax.set_ylim(lim)) for ax in pg.axes.flat[:]];
    # pg.fig.savefig('test.pdf')

    return ax


sort_order = \
['26N_pDNA',
 'D458_ETP',
 'D458_DMSO-1',
 'D458_DMSO-2',
 'D458_DMSO-3',
 'D458_DMSO-4',
 'D458_DMSO-5',
 'D458_JQ1-1',
 'D458_JQ1-2',
 'D458_JQ1-3',
 'D458_JQ1-4',
 'D458_JQ1-5',
 'PC9_ETP',
 'PC9_DMSO-1',
 'PC9_DMSO-2',
 'PC9_DMSO-3',
 'PC9_DMSO-4',
 'PC9_DMSO-5',
 'PC9_60nM-1',
 'PC9_60nM-2',
 'PC9_60nM-3',
 'PC9_60nM-4',
 'PC9_60nM-5',
 'PC9_1uM-1',
 'PC9_1uM-2',
 'PC9_1uM-3',
 'PC9_1uM-4',
 'PC9_1uM-5',
]