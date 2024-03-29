# (C) 2024, Tom Eulenfeld, MIT license

from sugar.data import gcode
import warnings


def biotranslate(seq, complete=False, gap='-', gap_after=2,
                 astop='X', warn=False, check_stop=False, tt=1):
    gc = gcode(tt)
    aas = []
    codon = ''
    ngap = 0
    check_start_codon = not complete
    for nb in str(seq).replace('U', 'T'):
        if nb == gap:
            ngap += 1
        else:
            codon = codon + nb
        if gap and gap_after is not None and ngap == gap_after:
            aas.append(gap)
            ngap -= 3
        if len(codon) == 3:
            if check_start_codon:
                if codon not in gc.starts and codon not in gc.astarts:
                    msg = (f'Codon {codon} is not a start codon {gc.starts} '
                           f'in genetic code #{gc.id}')
                    raise ValueError(msg)
                if warn and codon not in gc.starts:
                    msg = f'Codon {codon} possibly is not a start codon.'
                    warnings.warn(msg)
                check_start_codon = False
            try:
                aa = gc.tt[codon]
            except (KeyError):
                aa = 'X'
            if not complete and codon in gc.stops:
                break
            if codon in gc.astops:
                aa = astop
                if warn and not complete:
                    msg = f'Codon {codon} might be a stop codon.'
                    warnings.warn(msg)
            aas.append(aa)
            codon = ''
    else:
        if not complete and check_stop:
            msg = (f'Last codon is not a stop codon {gc.stops} '
                   f'in genetic code #{gc.id}')
            raise ValueError(msg)

    return ''.join(aas)
