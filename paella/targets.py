import regex as re
from collections import Counter, defaultdict
from itertools import product
from random import choice
import time
from paella.analysis import rc


random_spacer = 'ctactgtaccgatcaagcttc'
stop = 'TAG', 'TAA', 'TGA'

default_parts = {'sticky_5':  'AAAT',
         'sticky_3':  'GTAC',
         'PAM':       'cgg', # avoid start/stop codons
         'kozak_ATG': 'cgccaccatg',
         'BsmBI_site': 'CGTCTCc',
         'BsmBI_arm': 'tttcccggacgatc'
         }

sense = 'sticky_5', 'kozak_ATG', 'left_spacer', 'sgRNA', 'PAM', 'right_spacer'
antisense = 'sticky_5', 'kozak_ATG', 'left_spacer', 'rc_PAM', 'rc_sgRNA', 'right_spacer'

def build_oligo(design, custom_parts=None):
    parts = default_parts.copy()
    if custom_parts:
        parts.update(custom_parts)

    oligo = ''
    for element in design:
        if element.startswith('rc_'):
            oligo += rc(parts[element[3:]])
        else:
            oligo += parts[element]

    return oligo

def find_all(s, pat):
    return [m.start() for m in re.finditer(pat, s)]

def count_occurences(sequence, codon):
    a = [('sense', m.start() % 3) for m in re.finditer(codon, sequence)]
    b = [('antisense', m.start() % 3) for m in re.finditer(codon, rc(sequence))]
    return Counter(a + b)
    
def penalty_start(s):
    test = re.split('ATG', s.upper(), maxsplit=1)[1]
    penalty = len([x for x in find_all(test, 'ATG') if x % 3])
    return penalty

def penalty_stop(s):
    test = re.split('ATG', s.upper(), maxsplit=1)[1]
    penalty = sum(len(find_all(test, codon)) for codon in stop)
    return penalty

def score(s, activation_frame=1, cut_site_offset=0):
    """Count the number of 
    a) out-of-frame occurences of ATG in the sense strand, plus
    b) occurences of stop codons in activation frame after cut site.

    Automatically removes the first ATG and only analyzes subsequent 
    sequence.
    """
    test = re.split('ATG', s.upper(), maxsplit=1)[1]
    penalty = len([x for x in find_all(test, 'ATG') if x % 3])
    # penalty=0
    for codon in stop:
        # stop codon in the initial frame
        penalty += len([x for x in find_all(test, codon) if x % 3 == 0])
        # stop codon in the final frame
        penalty += len([x for x in find_all(test[cut_site_offset:], codon)
                        if x % 3 == activation_frame])

        return penalty

def get_frames(s, code):
    if code == 'start':
        xs = find_all(s, 'ATG')
    if code == 'stop':
        xs = [x for codon in stop for x in find_all(s, codon)]
    return sorted(set([x % 3 for x in xs]))
        
def get_targets(sg, spacer, translated_frames=(0, 1)):
    """Valid targets begin with a kozak and contain no ATG in the alternative frames.
    The target sequence (sgRNA and PAM) can be placed in the sense or antisense strand.
    """
    candidates = []
    for offset in range(len(spacer)):
        parts = {'sgRNA': sg, 'left_spacer': spacer[:offset], 'right_spacer': spacer[offset:]}
        PAM = default_parts['PAM']

        # sense
        candidate = build_oligo(sense, parts)
        cut_site_offset = candidate.find(sg + PAM) + len(sg + PAM) - 6 # - 3
        cut_site_offset -= 6
        # if score(candidate, cut_site_offset) == 0:
            # only perfect score
        if test_frames(candidate, translated_frames, stop_offset=cut_site_offset):
            candidates += [candidate]

        # antisense
        candidate = build_oligo(antisense, parts)
        cut_site_offset = candidate.find(rc(sg + PAM)) + 6 # - 3
        # if score(candidate, cut_site_offset) == 0:
            # only perfect score
        if test_frames(candidate, translated_frames, stop_offset=cut_site_offset):
            # print cut_site_offset, candidate
            candidates += [candidate]

    return candidates

def get_target_oligos(sgRNA, spacer):
    """Return top and bottom strand oligos ready for golden gate cloning.
    """
    sgRNA = sgRNA.upper()
    candidates = get_targets(sgRNA, spacer)
    arr = []
    for fwd in candidates:
        rev = rc(fwd[4:] + default_parts['sticky_3'])
        arr += [(fwd, rev)]
    return arr

def get_candidate_positions(sgRNA):
    targets = get_targets(sgRNA, random_spacer[:3])
    candidates = []
    for target in targets:
        start = target.upper().find('ATG')
        end = target.find(sgRNA)
        if end == -1:
            orientation = 'antisense'
            end = target.find(rc(sgRNA + default_parts['PAM']))
        else:
            orientation = 'sense'

        offset = end - start - 3
        candidates.append((orientation, offset)) 

    return candidates

def get_multiplexed_target(sgRNAs, translated_frames=(0, 1), attempts=1e4):
    """Target must contain no out-of-frame ATG and no stop codons in
    translated frames.
    """

    positions = [get_candidate_positions(s) for s in sgRNAs]

    bonus_parts = {'spacer_%d' % i: 'ccc'[:i] for i in range(3)}
    bonus_parts['random_spacer_front'] = random_spacer[:9]
    bonus_parts['random_spacer_back'] = random_spacer[-9:]

    positions_ = [[choice(x) for x in positions] for _ in range(int(attempts))]
    candidates = []

    for arr in positions_:
        target_arr = tuple()
        for sgRNA, (orientation, offset) in zip(sgRNAs, arr):
            bonus_parts['sgRNA_%s' % sgRNA] = sgRNA
            if orientation == 'sense': 
                insert = 'sgRNA_%s' % sgRNA, 'PAM'
            else:
                insert = 'rc_PAM', 'rc_sgRNA_%s' % sgRNA
            target_arr += 'spacer_%d' % offset,
            target_arr += insert
            target_arr += 'spacer_%d' % (3 - offset - 1),

        design = 'sticky_5', 'kozak_ATG', 'random_spacer_front'
        design += target_arr 
        design += 'random_spacer_back', 'sticky_3'
        extra_offset = (26 * len(arr)) % 3
        bonus_parts['random_spacer_back'] = random_spacer[-9:-extra_offset]
        candidate = build_oligo(design, bonus_parts)
        stops = get_frames(candidate, 'stop')
        if test_frames(candidate, translated_frames):
            candidates += [(score(candidate), candidate)]

    return sorted(candidates)


def test_frames(sequence, translated_frames, stop_offset=0):
    sequence = sequence.upper()
    start = get_frames(sequence, 'start')
    
    if len(start) > 1:
        return False
    else:
        start = start[0]
        # no stop codon in starting frame after ATG
        start_ix = sequence.find('ATG')
        stops = get_frames(sequence[start_ix:], 'stop')
        stops = [(s + start_ix) % 3 for s in stops]
        if start in stops:
            return False

        # no stop codon in other frames after stop offset
        stops = get_frames(sequence[stop_offset:], 'stop')
        stops = [(s + stop_offset) % 3 for s in stops]
        return not any((frame + start) % 3 in stops for frame in translated_frames)

def test_in_frame_stop(seq):    
    for codon in stop:
        for x in find_all(seq, codon):
            if x % 3 == 0:
                return True
    return False

def test_out_of_frame_start(seq):    
    for x in find_all(seq, 'ATG'):
        if x % 3 != 0:
            return True
    return False

def test_target(o_fwd, sgRNA):

    after_kozak = o_fwd.upper().split('GCCACCATG')[1]
    # has PAM
    if sgRNA in after_kozak:
        assert after_kozak.split(sgRNA)[1][1:3] == 'GG'
    else:
        assert after_kozak.split(rc(sgRNA))[0][-3:-1] == 'CC'

    # mimic 1 bp deletion
    if sgRNA in after_kozak:
        cut_site = after_kozak.index(sgRNA) + 17
    else:
        cut_site = after_kozak.index(rc(sgRNA)) + 3
        
    with_deletion = (after_kozak[:cut_site] + 
                     after_kozak[cut_site + 1:])

    a = test_in_frame_stop(with_deletion)
    b = test_in_frame_stop(after_kozak)
    c = test_out_of_frame_start(after_kozak)
    return not (a or b or c)