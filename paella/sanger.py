from paella.imports import *

import Bio.pairwise2
from Bio.pairwise2 import align


mNeon = """GTCAGCAAGGGCGAGGAGGACAACGCCTCTCTCCCAGCGA
CACACGAGTTACACATCTTTGGCTCCATCAACGGTGTGGACTT
CGACCTCGTGGGTCAGGGCACCGGCAATCCAAACGACGGTTAC
GAGGAGCTCAACCTCAAGTCCACCAAGGGCGACCTCCAGTTCT
CCCCCTGGATTCTGGTCCCTCATATCGGGTACGGCTTCCATCA
GTACCTGCCCTACCCAGACGGGATGTCGCCTTTCCAGGCCGCC
ATGGTGGACGGCTCCGGATACCAAGTCCATCGCACAATGCAGT
TCGAAGACGGTGCCTCCCTTACTGTCAACTACCGCTACACCTA
CGAGGGAAGCCACATCAAAGGAGAGGCCCAGGTCAAGGGGACT
GGTTTCCCTGCCGACGGTCCTGTCCTCACCAACTCGCTCACCG
CTGCGGACTGGTGCAGGTCGAAGAAGACTTACCCCAACGACAA
AACCATCATCAGTACCTTCAAGTGGAGTTACACCACTGGAAAC
GGCAAGCGCTACCGGAGCACTGCGCGGACCACCTACACCTTTG
CCAAGCCAATGGCGGCCAACTATCTCAAGAACCAGCCGATGTA
CGTGTTCCGCAAGACGGAGCTCAAGCACTCCAAGACCGAGCTC
AACTTCAAGGAGTGGCAAAAGGCCTTTACCGACGTCATGGGCA
TGGACGAGCTCTACAAG""".replace('\n', '')

def metrics(df_aligned):
    kozak = 'GCCACCATG'
    T2A = ('GGCAGTGGAGAGGGCAGAGGAAGTCT'
       'GCTCACCTGCGGCGACGTCGAGGAGAATCCTGGCCCA')
    k = 15
    T2A_kmers = [T2A[i:i+k] for i in range(len(T2A)-k + 1)]

    mNeon_kmers = [mNeon[i:i+k] for i in range(len(mNeon)-k + 1)]
    def score_mNeon(seq):
        for last_kmer in mNeon_kmers[::-1]:
            if last_kmer in seq:
                break
        else:
            return 0, 0
    
        kmers = mNeon_kmers[:mNeon_kmers.index(last_kmer)]
        score = np.mean([kmer in seq for kmer in kmers]) 
        coverage = len(kmers) / float(len(mNeon_kmers))
        return score, coverage
    
    def estimate_frameshift_short (seq, reporte):
        m = re.search('GCCACCATG(.+?)GGCAGTGGAGAGGGCAGAGGAAGTCTGCTCACCTGCGGCGACGTCGAGGAGAATCCTGGCCCA', reporter)
        if m:
            coding = m.group()
            coding_kmers = [coding[i:i+k] for i in range(len(coding)-k + 1)]

            try:
                seq = seq.split(kozak)[1]
            except IndexError:
                return -1, ''
            # only up until the first N
            seq = seq.split('N')[0]
            # find the last matching kmer
            for last_kmer_T2A in coding_kmers[::-1]:
                if last_kmer_T2A in seq:
                    break
            else:
                return -1, ''
            frameshift = seq.index(last_kmer_T2A) - coding.index(last_kmer_T2A)
            frameshift_T2A = frameshift % 3
        return frameshift_T2A, last_kmer_T2A
        
        
    def estimate_frameshift(seq, reporter): 
        coding = reporter.split(kozak)[1]
        coding_kmers = [coding[i:i+k] for i in range(len(coding)-k + 1)]
        
        # analyze read starting at kozak
        try:
            seq = seq.split(kozak)[1]
        except IndexError:
            return -1, ''
        # only up until the first N
        seq = seq.split('N')[0]
        # find the last matching kmer
        for last_kmer in coding_kmers[::-1]:
            if last_kmer in seq:
                break
        else:
            return -1, ''
        # frameshift between read and reference
        frameshift = seq.index(last_kmer) - coding.index(last_kmer)
        frameshift = frameshift % 3
        return frameshift, last_kmer

    def estimate_frameshift_T2A (seq, reporter):
        m = re.search('GCCACCATG(.+?)GGCAGTGGAGAGGGCAGAGGAAGTCTGCTCACCTGCGGCGACGTCGAGGAGAATCCTGGCCCA', reporter)
        if m:
            coding = m.group()
            coding_kmers = [coding[i:i+k] for i in range(len(coding)-k + 1)]

            try:
                seq = seq.split(kozak)[1]
            except IndexError:
                return -1, ''
                # only up until the first N
            seq = seq.split('N')[0]
            # find the last matching kmer
            for last_kmer_T2A in coding_kmers[::-1]:
                if last_kmer_T2A in seq:
                    break
            else:
                return -1, ''
            frameshift = seq.index(last_kmer_T2A) - coding.index(last_kmer_T2A)
            frameshift_T2A = frameshift % 3
            return frameshift_T2A, last_kmer_T2A    
    
    
    score_kmers = lambda x, kmers: sum(kmer in x for kmer in kmers)
    
    stop_3x = 'TAATAGTGAGTAGTAGTAAGTGATAATAG'
    arr = []
    for _, row in df_aligned.iterrows():
        mNeon_score, mNeon_coverage = score_mNeon(row['seq_al'])
        frameshift, last_kmer = estimate_frameshift(row['seq'], row['reporter_seq'])
        frameshift_T2A, last_kmer_T2A = estimate_frameshift_T2A (row['seq'], row['reporter_seq'])
        info = {'has_stop': stop_3x in row['seq_al'],
                'target_intact': row['target_seq'] in row['seq_al'],
                'T2A_score': score_kmers(row['seq_al'], T2A_kmers),
                'mNeon_score': mNeon_score,
                'mNeon_cov': mNeon_coverage,
                'read_length': approx_read_length(row['seq']),
                'matched_target': str(row['target']) in row['target_template'],
                'frameshift': frameshift,
                'last_kmer': last_kmer,
                'frameshift_T2A': frameshift_T2A,
                'last_kmer_T2A': last_kmer_T2A
               }
        arr += [info]
            
    return pd.DataFrame(arr)

def approx_read_length(seq):
    pat = '([^N]*)'
    try:
        return max(len(x) for x in re.findall(pat, seq))
    except ValueError:
        return 0


def load_seqs(files):
    pat = 'seq/(?P<target>...)[-_](?P<rep>\d+)[-_](?P<primer>.*).seq'
    df_sanger1 = (pd.DataFrame({'file': files})
     .assign(seq=lambda x: x['file'].apply(read_fasta))
     .pipe(lambda x: x.join(x['file'].str.extract(pat)))
     .dropna()
    )
    pat_2 = 'seq/(?P<target>...)[_][s][-_](?P<rep>\d+)[-_](?P<primer>.*).seq'
    df_sanger2 = (pd.DataFrame({'file': files})
     .assign(seq=lambda x: x['file'].apply(read_fasta))
     .pipe(lambda x: x.join(x['file'].str.extract(pat_2)))
     .dropna()
    )
    pat_3 = 'seq_2/(?P<target>...)[-_](?P<rep>\d+)[-_](?P<primer>.*).seq'
    df_sanger3 = (pd.DataFrame({'file': files})
     .assign(seq=lambda x: x['file'].apply(read_fasta))
     .pipe(lambda x: x.join(x['file'].str.extract(pat_3)))
    )
    
    df_sanger = pd.concat([df_sanger1, df_sanger2, df_sanger3])
    return df_sanger

def reporter_from_fasta(f):
    plasmid_seq = paella.analysis.read_fasta(f)[1]
    # trick for circular seq
    plasmid_seq = 3 * plasmid_seq.upper() 

    start = 'TGCAAAGATGGATAAAGTTT'
    # beginning of WPRE
    end = 'GTAATCAACCTCTGGATTAC'
    # beginning of blast
    end = 'GGATCCGGCGCAACAAACTTCTC'

    pat = '{0}(.*?){1}'.format(start, end)
    reporter = re.findall(pat, plasmid_seq.upper())[0]
    
    target_left = 'GATAATAGAAAT'
    target_right = 'GTACAGGCAGTGGA'

    pat = '{0}(.*?){1}'.format(target_left, target_right)
    target = re.findall(pat, reporter)[0]
    
    return {'file': f, 'reporter_seq': reporter, 'target_seq': target}

def align_seqs_targets(df_sanger, df_tm):
    arr = []
    it = list(product(df_sanger.iterrows(), df_tm.iterrows()))
    for (_, sanger_info), (_, target_info) in tqdn(it):
            
            alignment = align.globalms(sanger_info['seq'], target_info['reporter_seq'],
                                  2, -1, -.5, -.1)[0]

            seq_al, target_al, score, _, _ = alignment

            f = 'aligned_2/{template}_{target}_{rep}_{primer}'
            f = f.format(template=target_info['file'], **sanger_info.to_dict())
            txt = Bio.pairwise2.format_alignment(*alignment)
            with open(f, 'w') as fh:
                fh.write(txt)
                
            info = {'seq_al': seq_al, 'target_al': target_al, 'score': score,
                   'target_template': target_info['file'], 
                    'target_seq': target_info['target_seq'],
                   'reporter_seq': target_info['reporter_seq']}
            info.update(sanger_info.to_dict())
                
            arr += [info]
            
    return pd.DataFrame(arr)

