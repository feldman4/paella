from paella.imports import *

def plot_cdf(y, **kwargs):
    ax = plt.gca()
    x = np.linspace(0, 1, len(y))
    y = np.cumsum(y)
    y = y / y.max()
    ax.plot(x, y, **kwargs)
    
def plot_top_sg_clustermap(df2, df_info, top=100):
    x = (df2.sort_values('min_rank')
            .filter(regex='_frac', axis=1)
            .rename(columns=lambda x: x.replace('_frac', ''))
            .head(top)
            .pipe(np.log10))

    conditions = df_info.set_index('sample_info').loc[x.columns, 'condition']
    palette = sns.color_palette("Set2", len(set(conditions)))
    palette = {c: p for c,p in zip(set(conditions), palette)}
    row_colors = conditions.apply(lambda x: palette[x])

    cg = sns.clustermap(x.T, cbar_kws={'label': 'log10'}, row_colors=row_colors)
    cg.ax_heatmap.set_xticks([])
    cg.ax_heatmap.set_ylabel('')
    cg.ax_heatmap.set_xlabel('barcode')

    legend_TN = [matplotlib.patches.Patch(color=c, label=l) for l,c in palette.items()]
    # SO
    l2=cg.ax_heatmap.legend(loc='center left',bbox_to_anchor=(-0.35,0.85),handles=legend_TN,frameon=True)
    l2.set_title(title='condition',prop={'size':14})

    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    return cg

def plot_barcode_distributions(df):

    fg = \
    (df
        .pipe(sns.FacetGrid, hue='condition', col='complexity', sharex=False)
        .map(plt.plot, 'fraction')
        .add_legend())

    [ax.set_yscale('log') for ax in fg.axes.flatten()]

    spike_ins = df.groupby('complexity').head(1)['complexity_int']
    for ax, s in zip(fg.axes.flatten(), spike_ins):
        ax.set_xlim([0, s * 2])
        ax.set_xlabel('rank')

    plot_add_title(fg.fig, 'scaled by expected ETP barcodes')
    fg.savefig('20180123_AJG_DF/figures/distribution_1.pdf')
    display(fg.fig)

    spike_ins = df.groupby('complexity').head(1)['spike_in']
    for ax, s in zip(fg.axes.flatten(), spike_ins):
        ax.set_xlim([0, s * 2])
        ax.set_xlabel('rank')

    plot_add_title(fg.fig, 'scaled by expected hygro barcodes')
    fg.savefig('20180123_AJG_DF/figures/distribution_2.pdf')
    display(fg.fig)

    top_barcodes = 100
    [ax.set_xlim([0, top_barcodes]) for ax in fg.axes.flatten()];
    [ax.set_xlabel('rank') for ax in fg.axes.flatten()];
    plot_add_title(fg.fig, 'top barcodes')
    fg.savefig('20180123_AJG_DF/figures/distribution_3.pdf')

def plot_add_title(fig, title):
    fig.suptitle(title, 
                y=1.05, fontsize=14)
    fig.tight_layout()