from paella.imports import *


singlet_gate = and_join(['FSC_A > 2 & FSC_A < 6',
                 'FSC_H > 1.4 & FSC_H < 6'])

single_cell_gate = 'SSC_A > 4.5 & SSC_A < 6.5'

gates = {
    'singlets': singlet_gate,
    'single_cells': single_cell_gate
}

limits = {
    ECD: [0, 1],
    FITC: [0, 1]
}



def fix_row(df):
    """ duplicate row on plate1, Exp_20180322_1
    """
    filt = df['row'] == 'H'
    df.loc[filt, 'row'] = 'G'
    return df