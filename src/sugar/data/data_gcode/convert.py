# (C) 2024, Tom Eulenfeld, MIT license
import json
from itertools import product
from sugar.data import CODES


def generate_gc(id_, name, aas, special_codons):

    gc = {'id': id_, 'name': name, 'aa_line': aas, 'sc_line': special_codons}
    tt = {}
    starts = []
    stops = []
    for i, (b1, b2, b3) in enumerate(product('TCAG', repeat=3)):
        codon = b1+b2+b3
        tt[codon] = aas[i]
        if special_codons[i] == 'M':
            starts.append(codon)
        if special_codons[i] == '*':
            stops.append(codon)
    ttinv = {}
    for k, v in tt.items():
        ttinv.setdefault(v, []).append(k)
    # check for all codons which are possibly stop codons
    all_codes = set(CODES.keys()) - {'.', '-'}
    astops = []
    for i, (b1, b2, b3) in enumerate(product(all_codes, repeat=3)):
        codon = b1 + b2 + b3
        if codon in tt:
            continue
        if any(bb1+bb2+bb3  in stops for bb1 in CODES[b1]
               for bb2 in CODES[b2] for bb3 in CODES[b3]):
            astops.append(codon)
    # check for all codons which are possibly stop codons#
    astarts = []
    for i, (b1, b2, b3) in enumerate(product(all_codes, repeat=3)):
        codon = b1 + b2 + b3
        if codon in tt:
            continue
        if any(bb1+bb2+bb3 in starts for bb1 in CODES[b1]
               for bb2 in CODES[b2] for bb3 in CODES[b3]):
            astarts.append(codon)
    # check which codons with ambigous bases code for the same aa
    for i, (b1, b2, b3) in enumerate(product(all_codes, repeat=3)):
        codon = b1 + b2 + b3
        if codon in tt:
            continue
        if len(s := set(tt[bb1+bb2+bb3] for bb1 in CODES[b1]
                        for bb2 in CODES[b2] for bb3 in CODES[b3])
               ) == 1:
            tt[codon] = s.pop()
    gc.update({'tt': tt, 'ttinv': ttinv,
               'starts': starts, 'astarts': astarts,
               'stops': stops, 'astops': astops})
    return gc


def filter_line(line, s):
    return line.replace(s, '').replace(',', '').replace('"', '').strip()

gcs = {}

with open('gc.prt') as f:
    for line in f:
        if line.strip().startswith('name'):
            pname = filter_line(line, 'name')
            if not pname.startswith('SGC'):
                name = pname
        if line.strip().startswith('id'):
            id_ = int(filter_line(line, 'id'))
        elif line.strip().startswith('ncbieaa'):
            aas = filter_line(line, 'ncbieaa')
        elif line.strip().startswith('sncbieaa'):
            special_codons = filter_line(line, 'sncbieaa')
            gcs[id_] = generate_gc(id_, name, aas, special_codons)

with open('gc.json', 'w') as f:
    json.dump(gcs, f)
