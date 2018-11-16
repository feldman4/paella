from paella.imports import *
import paella.flow

singlet_gate = and_join(['FSC_A > 2 & FSC_A < 6',
                 'FSC_H > 1.4 & FSC_H < 6'])

single_cell_gate = 'SSC_A > 4.5 & SSC_A < 6.5'

gates = {
    'singlets': singlet_gate,
    'single_cells': single_cell_gate
}

limits = {
    ECD:   [3, 6],
    FITC:  [2.5, 6.5],
    APC:   [3, 6],
    PE:    [1.1, 5],
    FSC_A: [0, 6e5],
    SSC_A: [0, 1e6],
    FSC_H: [0, 3e5],
}


name_map = {'FL1_A': FITC,
 'FL10_A': PE,
 'FL3_A': APC,
 # 'inf': '9999999',
 }


PE_cutoff = 2.5
PE_factor = 4

def patch_PE_data(x):
    return paella.flow.compress_below(x, PE_cutoff, PE_factor)


def load_flowjo_gates(xml_file, gate_names):
    gates = paella.flow.xml_to_pandas_gates(xml_file, transform=None)
    gates = [fix_gate_names(g, name_map) for g in gates]
    flowjo_gates = {name: gate for name, gate in zip(gate_names, gates)}

    gate_tags = paella.flow.extract_polygon_gates(xml_file)
    coords = map(paella.flow.tag_to_coords, gate_tags)

    xy_coords = {name: np.array(list(xy) + [xy[0]]) for name, (xy, _) in zip(gate_names, coords)}

    return flowjo_gates, xy_coords


### flowjo gates

def fix_gate_names(gate, name_map):
    for a, b in name_map.items():
        gate = gate.replace(a,b)
    return gate


def fix_row(df):
    """ duplicate row on plate1, Exp_20180322_1
    """
    filt = df['row'] == 'H'
    df.loc[filt, 'row'] = 'G'
    return df

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

    # reporter_names = {'TM36': 'FR-EF1a-mNeon-Blast-mCherry', 
    #                 'TM42': 'FR-EF1a-Zeo-mNeon-Blast-mCherry', 
    #                 'TM43': 'FR-EF1a-H2K-mNeon-Blast-mCherr'
    #                  }

    return (df
         .assign(sgRNA=lambda x: x['sgRNA'].apply(sgRNA_names.get))      
         .assign(target=lambda x: x['target'].apply(target_names.get))
         # .assign(reporter=lambda x: x['reporter'].apply(reporter_names.get)) 
        )