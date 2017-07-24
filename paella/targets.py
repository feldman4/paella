
sticky_5 = 'AAAT'
kozac_ATG = 'cgccaccatg'
PAM = 'agg'
random_spacer = 'ctactgtagcgatcaagcttact'
sticky_3 = 'GTAC'

def form1(sg, offset): 
    return ''.join([sticky_5,
                    kozac_ATG,
                    sg,
                    random_spacer[:offset],
                    PAM,
                    random_spacer[offset:]])

def form2(sg, offset): 
    return ''.join([sticky_5,
                    kozac_ATG,
                    random_spacer[offset:],
                    rc(sg + PAM),
                    random_spacer[:offset]])

def find_all(s, pat):
    return [m.start() for m in re.finditer(pat, s)]

def optimize(sg):            
    def score(s):
        test = s.split(kozac_ATG)[1].upper()
        return len([x for x in find_all(test, 'ATG') if x%3])
        
    candidates = []
    for offset in range(3):
        for form in (form1, form2):
            candidates += [form(sg, offset)]
    
    candidates = [(score(c), c) for c in candidates]
    return sorted(candidates)

def targets(sgRNAs):
    arr = []
    for i, s in enumerate(sgRNAs):
        best_FWD = optimize(s)[0][1]
        best_REV = rc(best_FWD[4:] + sticky_3)
        arr += [[best_FWD, best_REV]]
    return arr
