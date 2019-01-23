from paella.imports import *
import holoviews as hv
import datashader as ds
import FlowCytometryTools
from matplotlib import path
from bs4 import BeautifulSoup
import paella


FSC_A_DIM = hv.Dimension(FSC_A, range=(3, 7))
FSC_H_DIM = hv.Dimension(FSC_H, range=(3.5, 5.8))
SSC_A_DIM = hv.Dimension(SSC_A, range=(2, 8))

FITC_DIM = hv.Dimension(FITC, range=(3.5, 6))
ECD_DIM  = hv.Dimension(ECD, range=(2, 6.3))

ROW_DIM = hv.Dimension('row', values=list('ABCDEFGH'))
COL_DIM = hv.Dimension('col', values=range(1, 13))
FCS_COLS = [FSC_A, SSC_A, FSC_H, FITC, ECD, TIME]


gates = {'mCherry'      : 'FITC_A < 4.5 & ECD_A > 4',
         'mNeon'        : 'FITC_A > 4.5 & ECD_A < 3.5',
         'double_neg'   : 'FITC_A < 4.5 & ECD_A < 3.5'
}


try:
    paella.gates.update(gates)
except AttributeError:
    paella.gates = gates

## specific to pipeline

def rename_variables_for_plot(df):
	sgRNA_names  = {'sg20N_304': 'guide 1',
	                'sg20N_305': 'guide X',
	                'sg20N_307': 'guide 2',
	                'sg20N_308': 'guide 3',
	                'sg20N_316': 'guide 4',
	                'sg20N_318': 'guide 5',
	                'sg20N_319': 'guide 6',
	                }

	target_names = {'T304': 'target 1', 
					'T305': 'target X', 
					'T307': 'target 2', 
					'T308': 'target 3', 
					'T316': 'target 4', 
					'T318': 'target 5',
					'T319': 'target 6',
					'MT0': 'Multi-target 1', 
					'MT1': 'Multi-target 2',
					}

	return (df
		 .assign(sgRNA=lambda x: x['sgRNA'].apply(sgRNA_names.get))      
		 .assign(target=lambda x: x['target'].apply(target_names.get))
		)

def plot_plate_summary(df, col_name, aggfunc='log_sum', **heatmap_kwargs):
    name = 'value'
    axs = []
    nrows = len(set(df['date']))
    ncols = len(set(df['plate']))
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*5, nrows*3), 
                            sharex=True, sharey=True)
    it = zip(axs.flatten(), df.groupby(['date', 'plate']))
        
    if aggfunc == 'mean':
        agg = lambda x: x.mean()
    elif aggfunc == 'log_sum':
        agg = lambda x: np.log10(1 + x.sum())
    elif aggfunc == 'log_mean':
        agg = lambda x: np.log10(1e-5 + x.mean())


    for ax, ((date, plate), df_) in it:
        z = \
        (df_
         .groupby(['row', 'col'])[col_name].pipe(agg)
         .rename(name).reset_index()
         .pivot_table(index='row', columns='col', values=name)
         .pipe(sns.heatmap, ax=ax, **heatmap_kwargs)
        )
        ax.set_title('{date} plate {plate}: {col_name}'.format(date=date, plate=plate, col_name=col_name))
        ax.set_xlabel('')
        ax.set_ylabel('')        
        
    return fig

def get_info(filename):
    pat = '(\d\d)-Tube-(.)(\d+).fcs'
    plate, row, col = re.findall(pat, filename)[0]
    return {'plate': int(plate), 'well': '%s%02d' % (row, int(col)), 
            'row': row, 'col': int(col)}

def load_file(filename, n):
    df = (FlowCytometryTools.FCMeasurement('', datafile=filename)
          .data
          [:n]
          .rename(columns=lambda x: x.replace('-', '_'))
          [FCS_COLS]
          .pipe(log_transform)
          .assign(file=filename)
          .assign(**get_info(filename))
         )
    
    return df

def log_transform(x):
    return np.log10(1 + np.abs(x))

def fix_row(df):
    """ duplicate row on plate1, Exp_20180322_1
    """
    filt = df['row'] == 'H'
    df.loc[filt, 'row'] = 'G'
    return df

def load_gate_files():
    gate_files = glob('TM_Figure2_FlowJo/*fcs_gates.xml')

    pat = 'Figure2_0(\d)-Tube-(.*).fcs_gates.xml'
    return (pd.DataFrame(map(lambda x: re.findall(pat, x)[0], gate_files))
     .rename(columns={0: 'plate', 1: 'well'})
     .assign(row=lambda x: x['well'].str.slice(stop=1))
     .assign(col=lambda x: x['well'].str.slice(start=1).astype(int))
     .assign(gate_file=gate_files)
    )

def make_layout():
    """Based on 20180323_Fig2_FACS/FlowLayout.xlsx
    """
    # reporter
    columns = ['TM36', 'TM42', 'TM43',
               'TM36', 'TM42', 'TM43',
               'TM36', 'TM42', 'TM43',
               'TM36', 'TM42', 'TM43']

    # sgRNA
    rows = ['sg20N_304', 'sg20N_305', 
            'sg20N_307', 'sg20N_308',
            'sg20N_316', 'sg20N_318',
            'sg20N_319']

    # target
    targets = ['T304', 'T305', 'T307', 'T308', 'T316', 'T318', 'T319']
    
    wells = ['%s%02d' % (r,c) for r in 'ABCDEFG' for c in range(1,13)]

    df = \
    (pd.DataFrame({'well': wells})
        .assign(row=lambda x: x['well'].apply(lambda y: y[0]))
        .assign(col=lambda x: x['well'].apply(lambda y: int(y[1:])))
        .assign(reporter=lambda x: x['col'].apply(lambda y: columns[y - 1]))
        .assign(sgRNA=lambda x: x['row'].apply(
            lambda y: rows['ABCDEFGH'.index(y)]))
        )

    def row_col_to_target(row, col):
        if col in (1, 2, 3):
            return targets['ABCDEFG'.index(row)]
        if col in (4, 5, 6):
            ix = ('ABCDEFG'.index(row) + 1) % 7
            return targets[ix]
        if col in (7, 8, 9):
            return 'MT0'
        if col in (10, 11, 12):
            return 'MT1'
        
    df['target'] = [row_col_to_target(r,c) for r,c in zip(df['row'], df['col'])]

    df1 = df.copy().assign(plate=1, dox=True)
    df2 = df.copy().assign(plate=2, dox=False)

    df = pd.concat([df1, df2]).query('sgRNA != "blank"')
    return df

def add_mNeon_mCherry_gate(df):
    flags = '1 * mNeon + 2 * mCherry + 4 * double_neg'
    flags = '1 * mNeon + 2 * mCherry' 
    categories = {0: 'double_pos', 1: 'mNeon', 
                  2: 'mCherry', 4: 'double_neg'}
    df['gate'] = (df.eval(flags).astype('category')
                   .cat.rename_categories(categories))
    return df

## holoviews

def apply_uni_rows(df):
    
    def uni_it(xs):
        return [unichr(60000 + i) + x for i,x in enumerate(xs)]

    rows = sorted(set(df['row']))[::-1]
    uni_map = OrderedDict(zip(rows, uni_it(rows)))
    df['row_uni'] = df['row'].apply(lambda x: uni_map[x])
    
    return df

def to_grid_space(df, kdims_inner, kdims_outer, element=hv.Points):

    gb = df.groupby(kdims_outer)
    x = [(k, element(v, kdims=kdims_inner)) for k,v in gb]   

    # GridSpace bug https://github.com/ioam/holoviews/issues/1787
    grid = hv.GridSpace(x, kdims=kdims_outer)

    return grid

def make_gate_callback(df, gate_name='gate', kdims=(FSC_A, SSC_A)):
    kdims = list(kdims)

    def callback(data):
        poly = data
        
        # debug
        hv.poly = poly
        
        try:   
            for i, (xs, ys) in enumerate(zip(poly['x'], poly['y'])):
          
                points = np.array(zip(xs, ys))

                p = path.Path(points)
                mask = p.contains_points(df[kdims])
                
                df[gate_name] = mask 
                df[gate_name] = df[gate_name].astype(bool)

                # stash in global dict
                gate = gates_from_poly(poly, kdims)
                try:                
                    paella.gates[gate_name] = gate
                except KeyError:
                    paella.gates = {gate_name: gate}
                break

        except Exception as e:
            hv.error = e
            
        return hv.Points(df, kdims=kdims)

    return callback

def raster_points(df, kdims, color_dim=None, normalization='log'):
    from holoviews.operation.datashader import datashade, aggregate, dynspread

    if color_dim:
        # copy because datashader.count_cat mutates datatype 
        # to categorical...
        columns = []
        for x in list(kdims) + [color_dim]:
            if isinstance(x, basestring):
                columns += [x]
            elif isinstance(x, hv.Dimension):
                columns += [x.label]
            else:
                raise ValueError
        points = hv.Points(df[columns].copy(), kdims)
        datashaded = datashade(points, 
                               aggregator=ds.count_cat(color_dim),
                               normalization=normalization)
    else:
        points = hv.Points(df, kdims)
        datashaded = datashade(points, 
                       normalization=normalization)
        
    return dynspread(datashaded)

